% THIS SCRIPT IS ONLY INTENDED FOR USE BY THE PACKAGE MAINTAINER

% Load required packages
pkg load generate_html
pkg load parallel
pkg load statistics
pkg load statistics-resampling

% Get options
myopt = get_html_options ("octave-forge");

% Customize the options
myopt.include_package_news = false;
myopt.include_alpha = false;
myopt.include_demos = false;
myopt.download_link = "https://gnu-octave.github.io/packages/statistics-resampling/";
myopt.repository_link = "https://github.com/gnu-octave/statistics-resampling/tree/master";
myopt.older_versions_download = "https://github.com/gnu-octave/statistics-resampling/releases";
generate_package_html ('statistics-resampling', 'manual_docs', myopt);

% Modify the footer
system ('sed -i '''' -e  ''s/\(.*\)index.html/\1overview.html/'' ./manual_docs/statistics-resampling/function/*.html');
system ('sed -i '''' -e  ''/footer.js/d'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''/footer.js/d'' ./manual_docs/statistics-resampling/function/*.html');

% Blank navigation pane
system ('sed -i '''' -e  ''/fixed.js/d'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''/fixed.js/d'' ./manual_docs/statistics-resampling/function/*.html');

% Delete files that are surplus to requirements
delete ('./manual_docs/footer.js')
delete ('./manual_docs/fixed.js')
delete ('./manual_docs/news.png')
delete ('./manual_docs/manual.png')

% Replace the javascript file
delete ('./manual_docs/javascript.js')
copyfile ('./javascript.js','./manual_docs/javascript.js')

% Create accessible link 
if exist ('./manual.html', 'file')
  delete ('./manual.html')
end
system ('ln -s ./manual_docs/statistics-resampling/index.html ./manual.html');

% Clear uplus
clear myopt
