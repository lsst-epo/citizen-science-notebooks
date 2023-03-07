---
name: Technical Issue
about: Report a technical issue/bug that you have observed while using the Citizen Science Notebooks in the RSP Notebook Aspect.
labels: bug
---

**Describe the bug**
A description of what the technical issue is. Some questions to consider addressing: How was the issue discovered? When did you first notice the bug? Is it affecting new functionality or existing, previously working functionality?

**To Reproduce**
Steps to reproduce the behavior, written in imperative mood:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Actual behavior**
A clear and concise description of what actually happened.

**Screenshots**
If applicable, add screenshots to help explain your problem. You do not need to post a screenshot of the output from the `send_data()` cell of the cSci notebook - please use the next section for that.

**EDC Output**
If applicable, copy and paste the output from the `Send the...to Zooniverse` cell into the below tilde wrapped box:

```
SAMPLE OUTPUT
'1. Checking batch status'
'    ** Warning! - The Zooniverse client is throwing an error about a missing subject set, this can likely safely be ignored.'
'2. Writing metadata file required by the Rubin EPO Data Center.'
'3. Zipping up all the astro cutouts - this can take a few minutes with large data sets, but unlikely more than 10 minutes.'
'4. Uploading the citizen science data'
'5. Creating a new Zooniverse subject set'
'6. Notifying the Rubin EPO Data Center of the new data, which will finish processing of the data and notify Zooniverse'
'** Additional information:'
root WARNING: Your project has not been approved by the data rights panel as of yet, as such you will not be able to send any additional data to Zooniverse until your project is approved.
'7. Sending the manifest URL to Zooniverse'
'** Information: subject_set.id: 111630; manifest: https://storage.googleapis.com/citizen-science-data-public/4de53abc-efbd-494c-ac13-9bfaff672329/manifest.csv'
'8. Transfer process complete, but further processing is required on the Zooniverse platform and you will receive an email at erosas@lsst.org'
```

**Additional context**
Add any other context about the problem here. 