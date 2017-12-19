mkdir -p $HOME/celgene-rnd-riku-researchanalytics
s3fs \
celgene-rnd-riku-researchanalytics $HOME/celgene-rnd-riku-researchanalytics \
-o url=https://s3.amazonaws.com \
-o use_sse=1 \
-o umask=0002 \
-o allow_other \
-o use_path_request_style \
-o uid=$(id -u) \
-o gid=$(id -g) \
-o passwd_file=~/.aws/passwd-s3fs 
