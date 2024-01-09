#!/bin/bash

# Delete files that are surplus to requirements
rm ./tmp/footer.js
rm ./tmp/news.png
rm ./tmp/manual.png
rm ./tmp/oct.png

# Move custom javascript and css documents to a subfolder in the dir containing index.html
mkdir ./tmp/statistics-resampling/site-files
cp ./templates/javascript.js ./tmp/statistics-resampling/site-files/javascript.js
cp ./templates/manual.css ./tmp/statistics-resampling/site-files/manual.css
cp ../doc/icon_48x48.png ./tmp/statistics-resampling/site-files/pkg_icon_48x48.png
rm ./tmp/javascript.js
rm ./tmp/octave-forge.css
mv ./tmp/fixed.js ./tmp/statistics-resampling/site-files/fixed.js
mv ./tmp/favicon.ico ./tmp/statistics-resampling/site-files/favicon.ico
mv ./tmp/doc.png ./tmp/statistics-resampling/site-files/doc.png
mv ./tmp/homepage.png ./tmp/statistics-resampling/site-files/homepage.png
mv ./tmp/repository.png ./tmp/statistics-resampling/site-files/repository.png
mv ./tmp/download.png ./tmp/statistics-resampling/site-files/download.png
mv ./tmp/statistics-resampling/COPYING.html ./tmp/statistics-resampling/copying.html
rm -Rf ./tmp/icons

# Update file references
sed -i '' -e  "s/write_top_menu (\(.*\));/write_top_menu ('\.');/g"  ./tmp/statistics-resampling/*.html
sed -i '' -e  "s/write_top_menu (\(.*\));/write_top_menu ('\.\.');/g" ./tmp/statistics-resampling/function/*.html
sed -i '' -e  "s/write_docs_left_menu (\(.*\));/write_docs_left_menu ('\.');/g" ./tmp/statistics-resampling/*.html
sed -i '' -e  "s/write_docs_left_menu (\(.*\));/write_docs_left_menu ('\.\.');/g" ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/octave-forge\.css/\.\/site-files\/manual\.css/g' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/\.\.\/octave-forge\.css/site-files\/manual\.css/g' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/fixed\.js/\.\/site-files\/fixed\.js/g' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/\.\.\/fixed\.js/site-files\/fixed\.js/g' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/javascript\.js/\.\/site-files\/javascript\.js/g' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/\.\.\/javascript\.js/site-files\/javascript\.js/g' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  '/footer\.js/d' ./tmp/statistics-resampling/*.html
sed -i '' -e  '/footer\.js/d' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/<script> write_footer (); <\/script>//g' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/<script> write_footer (); <\/script>//g' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/download\.png/\.\/site-files\/download\.png/g' ./tmp/statistics-resampling/index.html
sed -i '' -e  's/\.\.\/repository\.png/\.\/site-files\/repository\.png/g' ./tmp/statistics-resampling/index.html
sed -i '' -e  's/\.\.\/doc\.png/\.\/site-files\/doc\.png/g' ./tmp/statistics-resampling/index.html
sed -i '' -e  's/\.\.\/homepage\.png/\.\/site-files\/homepage\.png/g' ./tmp/statistics-resampling/index.html
sed -i '' -e  's/COPYING\.html/copying\.html/g' ./tmp/statistics-resampling/index.html
sed -i '' -e  's/\.\.\/favicon.ico"/\.\/site-files\/favicon.ico"/g' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/\.\.\/favicon.ico"/\.\/site-files\/favicon.ico"/g' ./tmp/statistics-resampling/function/*.html

# Change page headings of the index (package info) and function reference pages etc.
sed -i '' -e  's/statistics-resampling<\/h2>/About this package<\/h2>/g' ./tmp/statistics-resampling/index.html
sed -i '' -e  's/Octave<\/a> >= 4\.4\.0/Octave<\/a> >= 4\.4\.0 or Matlab >= R2007a 7\.4\.0/g' ./tmp/statistics-resampling/index.html
sed -i '' -e  's/statistics-resampling<\/h2>/Function Reference<\/h2>/g' ./tmp/statistics-resampling/function_reference.html
sed -i '' -e  's/<h2>/<h3>/g' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/<\/h2>/<\/h3>/g' ./tmp/statistics-resampling/function/*.html
cd ./tmp/statistics-resampling/function/
for f in *.html; do sed -i '' -e "s/<pre>/<h2>$(echo "${f%.*}")<\/h2>\n<pre>/" $f; done
for i in {1..20}; do sed -i '' -e 's/<h3>Demonstration '$i'/\n<h3><a name="'$i'">Demonstration '$i'<\/a>/' ./tmp/statistics-resampling/function/*.html; done
cd ../../..

# Publish README markdown file as HTML and add header
cd ..
pandoc README.md -o ./man/tmp/statistics-resampling/tmp.html
cat ./man/templates/header.html ./man/tmp/statistics-resampling/tmp.html > ./man/tmp/statistics-resampling/readme.html
printf "</div>\n</body>\n</html>\n" >> ./man/tmp/statistics-resampling/readme.html
rm ./man/tmp/statistics-resampling/tmp.html
cd man

# Add hyperlink to function reference directly in the HTML of the readme file (since the link in the markdown file would fail on the GitHub page) 
#sed -i '' -e  's/Function Reference/\<a href\=\"function_reference\.html\"\>Function Reference\<\/a\>/g' ./tmp/statistics-resampling/readme.html

# Update meta information in all html files
sed -i '' -e  's/<meta name="author" content="\(.*\)" \/>/<meta name="author" content="Andrew Penn" \/>/' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/<meta name="author" content="\(.*\)" \/>/<meta name="author" content="Andrew Penn" \/>/' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/<title>\(.*\)<\/title>/<title>The statistics-resampling package manual<\/title>/' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/<title>\(.*\)<\/title>/<title>The statistics-resampling package manual<\/title>/' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/<meta name="description" content="\(.*\)" \/>/<meta name="description" content="A package for statistical analysis using resampling methods." \/>/' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/<meta name="description" content="\(.*\)" \/>/<meta name="description" content="A package for statistical analysis using resampling methods." \/>/' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/<meta name="keywords" lang="en" content="\(.*\)" \/>/<meta name="keywords" lang="en" content="GNU Octave Packages, MATLAB Toolbox\" \/>/' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/<meta name="keywords" lang="en" content="\(.*\)" \/>/<meta name="keywords" lang="en" content="GNU Octave Packages, MATLAB Toolbox\" \/>/' ./tmp/statistics-resampling/function/*.html

# Fixes for validation of XHTML 1.0 Strict
sed -i '' -e  's/<script>/<script type="text\/javascript">/g' ./tmp/statistics-resampling/*.html
sed -i '' -e  's/<script>/<script type="text\/javascript">/g' ./tmp/statistics-resampling/function/*.html
sed -i '' -e  's/alt="Repository icon">/alt="Repository icon"\/>/' ./tmp/statistics-resampling/index.html

# Copy web pages and files to docs folder and clean up
cp -R ./tmp/statistics-resampling/* ../docs/
rm -Rf ./tmp
