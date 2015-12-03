all: html pdf

pdf:
	pandoc -S writeup.md -o writeup.pdf

html:
	pandoc writeup.md --to html5 -s --filter pandoc-crossref --filter pandoc-citeproc -o writeup.html --mathjax

watch:
	reload -b -s writeup.html & make html & fswatch -o writeup.md | xargs -n1 -I% make html
