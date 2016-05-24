egrep $1 $2 | awk '{if (/^\w/){ print ">" $1 "\n" $10 "\n"}}'
