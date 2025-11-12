# Privacy and Data Handling

This Shiny application processes user-provided files only within the active session.  
All uploads are handled in temporary storage and automatically deleted after the session ends.  
No uploaded files are retained, logged, or transmitted outside the container environment.

### Deployment context
The app is hosted under **https://r-odaf.nl** inside a secure Docker container managed by the R-ODAF team at Maastricht University.  
Only standard server access logs (IP, timestamp, browser type) are stored briefly for technical monitoring and rotated regularly.

### User responsibility
Users should avoid uploading personally identifiable or confidential data.  
This instance is designed solely for research-grade RNA-seq count matrices and metadata files.
