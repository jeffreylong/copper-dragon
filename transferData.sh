#!/bin/bash

# script to transfer data from one s3 location to another
# filename1 is the name of the file to transfer (with its full s3 location)
# destination is the 'directory' on s3 where the file will be copied
# paired can be yes or no.
#               If yes then a second file will be assumed with a filename like the first but different mate identification.
#               If no then only that file is copied

if [ -z "$1" ] ; then
	echo "$0 <object list to transfer> <destination s3> <paired (yes|no)>"
	exit 0
fi
suffix='_merged.fastq'
input=$1
destination="$2"
paired=$3
if [ "$paired" == "no"  ] ; then
    paired="no"
elif [ "$paired" == "yes" ] ; then
    paired="yes"
else
    echo "Cannot recognize the paired option (third argument)"
    exit
fi



while IFS= read -r filename1
do
	echo "Preparing a script for data transfer of file \"$filename1\""
	source=$(dirname "$filename1")
	if [ "$paired" == "yes" ] ;then filename2=$(echo "$filename1" | sed "s/1${suffix}/2${suffix}/"); fi
	# get a random prefix
	rand=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
	
	sfn=$(basename "$filename1").bsub
	echo "#!/bin/bash" > $sfn
	echo "#BSUB -q spot" >> $sfn
	
	
	echo "Filename 1: \"$filename1\""
	echo "Filename 2: \"$filename1\""
	
	
	for fn in "$filename1" "$filename2" ; do
	    if [ -z "$fn" ] ; then break; fi
	    echo "Now setting up $fn"
		# get the data
	    echo "#Getting the file $fn from S3"  >> $sfn
	    echo "filename=\$(basename \"$fn\")" >> $sfn
	    echo "aws s3 cp --profile datadrop \"$source/\$filename\" \"\$NGS_TMP_DIR/\$filename\"" >> $sfn
	
		# upload the file
	    echo "#Uploading file $filename to S3" >> $sfn
	    echo "aws s3 cp --sse AES256 --profile ra \"\$NGS_TMP_DIR/\$filename\" \"$destination$rand-\$filename\"" >> $sfn
		echo "rm \"\$NGS_TMP_DIR/\$filename\"" >> $sfn
		# cleanup
	
	done

done < $input
