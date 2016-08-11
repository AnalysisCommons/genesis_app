#!/bin/bash 
md5_checksum_file="$1"
archive_directory="$2"
dna_nexus_archives="Commons:/Tools/Archive/"


if [ ! -e ${md5_checksum_file} ]
then
   echo "File  ${md5_checksum_file} is missing"
   exit 1
fi

if [ ! -d ${archive_directory} ]
then
   echo "Directory  ${archive_directory} is missing"
   exit 1
fi



if md5sum -c --status ${md5_checksum_file}
then
    echo "Archives don't need to be updated"
    exit 0
else
    archive_files_to_update=$( basename -a $(md5sum -c ${md5_checksum_file} | grep FAILED | awk -F: '{print $1}'))
fi



if dx --help  &> /dev/null;
then
    if dx ls "${dna_nexus_archives}" &> /dev/null;
    then
        for archive_file in ${archive_files_to_update}
        do
            dx download "${dna_nexus_archives}/${archive_file}" -o "${archive_directory}"
        done
    else
      echo "No access to the DNA NEXUS ARCHIVE DIR: ${dna_nexus_archives}"
      exit 1 
    fi
else
   echo "Install DNA NEXUS DX toolkit and rerun"
   exit 1
fi

