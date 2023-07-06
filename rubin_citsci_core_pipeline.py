# Author: Eric Rosas
# Note: DO NOT MODIFY!
# To-do: Find a "do not modify" and/or "no warranty use as-is" license to insert here

import csv, uuid, os, shutil, json, logging, urllib.request
from datetime import datetime, timezone, timedelta
from IPython.display import display
from google.cloud import storage
import panoptes_client
from panoptes_client import Project, SubjectSet, Classification

class CitSciPipeline:
    
    def __init__(self):
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/opt/lsst/software/jupyterlab/butler-secret/butler-gcs-idf-creds.json"
        self.vendor_batch_id = 0
        self.project_id = -1
        self.guid = ""
        self.manifest_url = ""
        self.edc_response = ""
        self.step = 0
        self.email = ""
        self.project = None
        self.client = None

    def login_to_zooniverse(self, slug_name, email):
        self.client = panoptes_client.Panoptes.connect(login="interactive")
        self.project = Project.find(slug=slug_name)
        self.project_id = self.project.id
        self.email = email
        print("You now are logged in to the Zooniverse platform.")
        return

    def write_metadata_file(self, manifest, batch_dir):    
        manifest_filename = 'manifest.csv'
        with open(batch_dir + manifest_filename, 'w', newline='') as csvfile:
            fieldnames = list(manifest[0].keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for cutout in manifest:
                writer.writerow(cutout)

        return f"{batch_dir}{manifest_filename}"
    
    def clean_up_unused_subject_set(self):
        self.log_step("Cleaning up unused subject set on the Zooniverse platform, vendor_batch_id : " + str(self.vendor_batch_id))

        try:
            subject_set = SubjectSet.find(str(self.vendor_batch_id))

            if subject_set.id == self.vendor_batch_id:
                subject_set.delete()

        except:
            display(f"** Warning: Failed to find the subject set with id: {str(self.vendor_batch_id)}- perhaps it's been deleted?.")
        return

    def send_zooniverse_manifest(self):
        self.log_step("Sending the manifest URL to Zooniverse")
        display("** Information: subject_set.id: " + str(self.vendor_batch_id) + "; manifest: " + self.manifest_url);

        payload = {"subject_set_imports": {"source_url": self.manifest_url, "links": {"subject_set": str(self.vendor_batch_id)}}}
        json_response, etag = self.client.post(path='/subject_set_imports', json=payload)
        return

    def create_new_subject_set(self, name):
        self.log_step("Creating a new Zooniverse subject set")

        # Create a new subject set
        subject_set = panoptes_client.SubjectSet()
        subject_set.links.project = self.project

        # Give the subject set a display name (that will only be visible to you on the Zooniverse platform)
        subject_set.display_name = name 
        subject_set.save()
        self.project.reload()
        self.vendor_batch_id = subject_set.id
        return self.vendor_batch_id

    def check_status(self):
        status_uri = "https://rsp-data-exporter-dot-skyviewer.uw.r.appspot.com/citizen-science-ingest-status?guid=" + self.guid
        raw_response = urllib.request.urlopen(status_uri).read()
        response = raw_response.decode('UTF-8')
        return json.loads(response)

    def download_batch_metadata(self):
        project_id_str = str(self.project_id)
        dl_response = "https://rsp-data-exporter-dot-skyviewer.uw.r.appspot.com/active-batch-metadata?vendor_project_id=" + project_id_str
        raw_response = urllib.request.urlopen(dl_response).read()
        response = raw_response.decode('UTF-8')
        return json.loads(response)


    # Validates that the RSP user is allowed to create a new subject set
    def send_image_data(self, subject_set_name, batch_dir, cutout_data = None):
        self.step = 0
        self.log_step("Checking batch status")
        if self.has_active_batch() == True:
            raise CitizenScienceError("You cannot send another batch of data while a subject set is still active on the Zooniverse platform - you can only send a new batch of data if all subject sets associated to a project have been completed.")
        zip_path = self.zip_hips_cutouts(batch_dir)
        self.upload_cutouts(zip_path)
        self.create_new_subject_set(subject_set_name)

        self.edc_response = self.alert_edc_of_new_citsci_data()
        if(self.edc_response == None):
            self.edc_response = { "status": "error", "messages": "An error occurred while processing the data transfer process upload" }
        else:
            self.edc_response = json.loads(self.edc_response)

        if self.edc_response["status"] == "success":
            self.manifest_url = self.edc_response["manifest_url"]
            if len(self.edc_response["messages"]) > 0:
                display("** Additional information:")
                for message in self.edc_response["messages"]:
                    logging.warning(message)
                    # display("    ** " + message)
            else:
                self.log_step("Success! The URL to the manifest file can be found here:")
                display(self.manifest_url)
        else:
            self.clean_up_unused_subject_set()
            logging.error("** One or more errors occurred during the last step **")
            logging.error(self.edc_response["messages"])
            logging.error(f"Email address: {self.email}")
            logging.error(f"Timestamp: {str(datetime.now(timezone(-timedelta(hours=7))))}")
            # for message in edc_response["messages"]:
            #     display("        ** " + message)
            return

        self.send_zooniverse_manifest()
        self.log_step("Transfer process complete, but further processing is required on the Zooniverse platform and you will receive an email at " + self.email)
        return
    
    def zip_hips_cutouts(self, batch_dir):
        self.guid = str(uuid.uuid4())
        self.log_step("Zipping up all the astro cutouts - this can take a few minutes with large data sets, but unlikely more than 10 minutes.")
        shutil.make_archive("./" + self.guid, 'zip', batch_dir)
        return ["./" + self.guid + '.zip', self.guid + '.zip']

    def upload_cutouts(self, zip_path):
        self.log_step("Uploading the citizen science data")
        bucket_name = "citizen-science-data"
        destination_blob_name = zip_path[1]
        source_file_name = zip_path[0]

        storage_client = storage.Client()
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(destination_blob_name)

        blob.upload_from_filename(source_file_name)
        return

    def alert_edc_of_new_citsci_data(self):
        project_id_str = str(self.project_id)
        self.log_step("Notifying the Rubin EPO Data Center of the new data, which will finish processing of the data and notify Zooniverse")

        try:
            edc_endpoint = "https://rsp-data-exporter-dot-skyviewer.uw.r.appspot.com/citizen-science-bucket-ingest?email=" + self.email + "&vendor_project_id=" + project_id_str + "&guid=" + self.guid + "&vendor_batch_id=" + str(self.vendor_batch_id) + "&debug=True"
            response = urllib.request.urlopen(edc_endpoint).read()
            manifestUrl = response.decode('UTF-8')
            return manifestUrl
        except Exception as e:
            self.clean_up_unused_subject_set()
            return None

    # def send_butler_data_to_edc():
    #     log_step("Notifying the Rubin EPO Data Center of the new data, which will finish processing of the data and notify Zooniverse")
    #     edcEndpoint = "https://rsp-data-exporter-e3g4rcii3q-uc.a.run.app/citizen-science-butler-ingest?email=" + email + "&collection=" + datasetId + "&sourceId=" + sourceId + "&vendorProjectId=" + str(projectId) + "&vendor_batch_id=" + str(vendor_batch_id)
    #     log_step('Processing data for Zooniverse, this may take up to a few minutes.')
    #     response = urllib.request.urlopen(edcEndpoint).read()
    #     manifestUrl = response.decode('UTF-8')
    #     return

    def has_active_batch(self):
        active_batch = False
        for subject_set in self.project.links.subject_sets:
            try:
                for completeness_percent in list(subject_set.completeness.values()):
                    if completeness_percent == 1.0:
                        active_batch = True
                        break
                if active_batch:
                    break
            except:
                display("    ** Warning! - The Zooniverse client is throwing an error about a missing subject set, this can likely safely be ignored.");
        return active_batch

    def log_step(self, msg):
        self.step += 1
        display(str(self.step) + ". " + msg)
        return

    # Custom error handling for this notebook
    class CitizenScienceError(Exception):

        # Constructor or Initializer
        def __init__(self, value):
            self.value = value

        # __str__ is to print() the value
        def __str__(self):
            return(repr(self.value))
