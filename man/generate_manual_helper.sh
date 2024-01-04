#!/bin/bash

# Delete files that are surplus to requirements
rm ./manual_docs/footer.js
rm ./manual_docs/news.png
rm ./manual_docs/manual.png
rm ./manual_docs/oct.png
rm ./manual_docs/favicon.ico

# Move custom javascript and css documents to a subfolder in the dir containing index.html
mkdir ./manual_docs/statistics-resampling/site-files
cp ./templates/javascript.js ./manual_docs/statistics-resampling/site-files/javascript.js
cp ./templates/manual.css ./manual_docs/statistics-resampling/site-files/manual.css
cp ../doc/icon_48x48.png ./manual_docs/statistics-resampling/site-files/pkg_icon_48x48.png
rm ./manual_docs/javascript.js
rm ./manual_docs/octave-forge.css
mv ./manual_docs/fixed.js ./manual_docs/statistics-resampling/site-files/fixed.js
mv ./manual_docs/doc.png ./manual_docs/statistics-resampling/site-files/doc.png
mv ./manual_docs/homepage.png ./manual_docs/statistics-resampling/site-files/homepage.png
mv ./manual_docs/repository.png ./manual_docs/statistics-resampling/site-files/repository.png
mv ./manual_docs/download.png ./manual_docs/statistics-resampling/site-files/download.png
rm -Rf ./manual_docs/icons

# Update file references
sed -i '' -e  "s/write\_top\_menu \(\(.*\)\)\;/write\_top\_menu \(\'\.\'\)\;/g" ./manual_docs/statistics-resampling/*.html
sed -i '' -e  "s/write\_top\_menu \(\(.*\)\)\;/write\_top\_menu \(\'\.\.\'\)\;/g" ./manual_docs/statistics-resampling/function/*.html
sed -i '' -e  "s/write\_docs\_left\_menu \(\(.*\)\)\;/write\_docs\_left\_menu \(\'\.\'\)\;/g" ./manual_docs/statistics-resampling/*.html
sed -i '' -e  "s/write\_docs\_left\_menu \(\(.*\)\)\;/write\_docs\_left\_menu \(\'\.\.\'\)\;/g" ./manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/octave\-forge\.css/\.\/site\-files\/manual\.css/g' ./manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\.\.\/octave\-forge\.css/site\-files\/manual\.css/g' ./manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/fixed\.js/\.\/site\-files\/fixed\.js/g' ./manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\.\.\/fixed\.js/site\-files\/fixed\.js/g' ./manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/javascript\.js/\.\/site\-files\/javascript\.js/g' ./manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\.\.\/javascript\.js/site\-files\/javascript\.js/g' ./manual_docs/statistics-resampling/function/*.html
sed -i '' -e  '/footer\.js/d' ./manual_docs/statistics-resampling/*.html
sed -i '' -e  '/footer\.js/d' ./manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\<script\> write\_footer ()\; \<\/script\>//g' ./manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\<script\> write\_footer ()\; \<\/script\>//g' ./manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\.\.\/download\.png/\.\/site\-files\/download\.png/g' ./manual_docs/statistics-resampling/index.html
sed -i '' -e  's/\.\.\/repository\.png/\.\/site\-files\/repository\.png/g' ./manual_docs/statistics-resampling/index.html
sed -i '' -e  's/\.\.\/doc\.png/\.\/site\-files\/doc\.png/g' ./manual_docs/statistics-resampling/index.html
sed -i '' -e  's/\.\.\/homepage\.png/\.\/site\-files\/homepage\.png/g' ./manual_docs/statistics-resampling/index.html

# Change page headings of the index (package info) and function reference pages
sed -i '' -e  's/statistics\-resampling\<\/h2\>/About this package\<\/h2\>/g' ./manual_docs/statistics-resampling/index.html
sed -i '' -e  's/statistics\-resampling\<\/h2\>/Function Reference\<\/h2\>/g' ./manual_docs/statistics-resampling/function_reference.html

# Change to main package directory
cd ..

# Publish README markdown file as HTML and add header
pandoc README.md -o ./man/manual_docs/statistics-resampling/tmp.html
cat ./man/templates/header.html ./man/manual_docs/statistics-resampling/tmp.html > ./man/manual_docs/statistics-resampling/readme.html
printf "</div>\n</body>\n</html>\n" >> ./man/manual_docs/statistics-resampling/readme.html
rm ./man/manual_docs/statistics-resampling/tmp.html

# Add hyperlink to function reference directly in the HTML of the readme file (since the link in the markdown file would fail on the GitHub page) 
#sed -i '' -e  's/Function Reference/\<a href\=\"function_reference\.html\"\>Function Reference\<\/a\>/g' ./man/manual_docs/statistics-resampling/readme.html

# Update meta information in all html files
sed -i '' -e  's/\<meta name\=\"author\" content\=\"\(.*\)\" \/\>/\<meta name\=\"author\" content\=\"Andrew Penn\" \/\>/' ./man/manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\<meta name\=\"author\" content\=\"\(.*\)\" \/\>/\<meta name\=\"author\" content\=\"Andrew Penn\" \/\>/' ./man/manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\<title\>\(.*\)\<\/title\>/<title\>The statistics\-resampling package manual\<\/title\>/' ./man/manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\<title\>\(.*\)\<\/title\>/<title\>The statistics\-resampling package manual\<\/title\>/' ./man/manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\<meta name\=\"description\" content\=\"\(.*\)\" \/\>/\<meta name\=\"description\" content\=\"A package for statistical analysis using resampling methods.\" \/\>/' ./man/manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\<meta name\=\"description\" content\=\"\(.*\)\" \/\>/\<meta name\=\"description\" content\=\"A package for statistical analysis using resampling methods.\" \/\>/' ./man/manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/\<meta name\=\"keywords\" lang\=\"en\" content\=\"\(.*\)\" \/\>/\<meta name\=\"keywords\" lang\=\"en\" content\=\"GNU Octave Packages, MATLAB Toolbox\" \/\>/' ./man/manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\<meta name\=\"keywords\" lang\=\"en\" content\=\"\(.*\)\" \/\>/\<meta name\=\"keywords\" lang\=\"en\" content\=\"GNU Octave Packages, MATLAB Toolbox\" \/\>/' ./man/manual_docs/statistics-resampling/function/*.html
sed -i '' -e  '/\<link rel\=\"shortcut icon\" href\=\"\.\.\/favicon.ico\" \/\>/d' ./man/manual_docs/statistics-resampling/*.html
sed -i '' -e  '/\<link rel\=\"shortcut icon\" href\=\"\.\.\/.\.\/favicon.ico\" \/\>/d' ./man/manual_docs/statistics-resampling/function/*.html

# Fixes for validation of XHTML 1.0 Strict
sed -i '' -e  's/\<script\>/\<script type\=\"text\/javascript\"\>/g' ./man/manual_docs/statistics-resampling/*.html
sed -i '' -e  's/\<script\>/\<script type\=\"text\/javascript\"\>/g' ./man/manual_docs/statistics-resampling/function/*.html
sed -i '' -e  's/alt\=\"Repository icon\"\>/alt\=\"Repository icon\"\/\>/' ./man/manual_docs/statistics-resampling/index.html

# Create accessible link to in manual HTML pages
#rm ./manual
#ln -s ./man/manual_docs/statistics-resampling/ manual
mv ./man/manual_docs/statistics-resampling/* ./docs/
rm -Rf ./man/manual_docs

# Change back to man directory
cd man