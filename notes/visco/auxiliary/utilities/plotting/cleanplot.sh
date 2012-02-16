#!/bin/sh

perl -pi -e 's/BL { stroke userlinewidth 2/BL { stroke userlinewidth 0.2/' $1
perl -pi -e 's/LT0 { PL \[\] 1 0 0/LT0 { PL \[\] 0.76 0.08 0.02/' $1
perl -pi -e 's/LT1 { PL \[4 dl 2 dl\] 0 1 0/LT1 { PL \[4 dl 2 dl\] 0.00 0.33 0.54/' $1