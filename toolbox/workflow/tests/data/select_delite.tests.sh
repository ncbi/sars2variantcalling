#!/bin/bash -eu

( cat select_delite_tests.sql ; echo -e ";\ngo\n" ) | /opt/sybase/utils/bin/sqsh-ms -h -U anyone -P allowed -S SRA_PRIMARY1 -D SRA_Main > delite_tests.csv
