#!/bin/bash



gzip -d webdocs_1.data.gz
gzip -d webdocs_2.data.gz
gzip -d webdocs_3.data.gz
gzip -d webdocs_4.data.gz
gzip -d webdocs_5.data.gz
gzip -d webdocs_6.data.gz
gzip -d webdocs_7.data.gz

cat webdocs_1.data webdocs_2.data webdocs_3.data webdocs_4.data webdocs_5.data webdocs_6.data webdocs_7.data > webdocs.data




exit
