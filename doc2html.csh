#!/usr/bin/csh
if ($#argv < 1) then
	echo translate a multalin doc file to an Html file
	echo syntax: \'$0 myfile\' or \'$0 - \<myfile\'
	exit
endif
set separator = `grep -n '^//$' $1 |cut -f1 -d':'`
cat head.html
head -$separator $1
cat middle.html
@ separator++
more +$separator $1|sed '1,$s/)(//g;1,$s/\]\[//g;1,$s/(/<em class=low>/g;1,$s/\[/<em class=high>/g;1,$s/)/<\/em>/g;1,$s/\]/<\/em>/g'
cat tail.html
