#!/usr/bin/perl -w
die "translate a multalin doc file to an Html file\nsyntax: '$0 myinputfile' or '$0 - <myinputfile'\n" unless (@ARGV);
open (IN,$ARGV[0]) or die "Sorry! Cannot find the data\n";
print STDOUT <<END
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
 "http://www.w3.org/TR/html401">
<html>
<head>
<title>MultAlin result </title>
<style type='text/css'>
pre.seq {color: black; background-color: white}
em {font-style: normal}
em.low {color: blue; background-color: white}
em.high {color: red; background-color: white}
</style>
</head>
<body>
<pre><A HREF="#Alignment">Go directly to Alignment </A>
END
;
print STDOUT $_ while (defined($_=<IN>) && !/^\/\//);
print STDOUT "$_</pre><pre class=seq><A NAME='Alignment'></A>" ;
while (defined($_=<IN>))
{
        s/\)\(//g;
        s/\]\[//g;
        s/\(/<em class=low>/g;
        s/\[/<em class=high>/g;
        s/\)/<\/em>/g;
        s/\]/<\/em>/g;
        print STDOUT $_;
}
print STDOUT <<END
</pre>
</body>
</html>
END
;
