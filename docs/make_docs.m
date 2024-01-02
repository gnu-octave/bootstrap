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
myopt.include_demos = true;
myopt.download_link = "https://gnu-octave.github.io/packages/statistics-resampling/";
myopt.repository_link = "https://github.com/gnu-octave/statistics-resampling/tree/master";
myopt.older_versions_download = "https://github.com/gnu-octave/statistics-resampling/releases";
generate_package_html ('statistics-resampling', 'manual_docs', myopt);

% Modify the footer
system ('sed -i '''' -e  ''s/\(.*\)index.html/\1overview.html/'' ./manual_docs/statistics-resampling/function/*.html');
system ('sed -i '''' -e  ''s/\(.*\)index.html/\1overview.html/'' ./manual_docs/statistics-resampling/overview.html');
system ('sed -i '''' -e  ''s/\(.*\)statistics-bootstrap/\1Function Reference.html/'' ./manual_docs/statistics-resampling/function/*.html');
system ('sed -i '''' -e  ''s/\(.*\)statistics-bootstrap/\1Function Reference.html/'' ./manual_docs/statistics-resampling/overview.html');
system ('sed -i '''' -e  ''/footer.js/d'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''/footer.js/d'' ./manual_docs/statistics-resampling/function/*.html');

% Delete files that are surplus to requirements
delete ('./manual_docs/footer.js')
delete ('./manual_docs/news.png')
delete ('./manual_docs/manual.png')
delete ('./manual_docs/oct.png')

% Replace the javascript file
delete ('./manual_docs/javascript.js')
delete ('./manual_docs/octave-forge.css')
copyfile ('./templates/javascript.js','./manual_docs/javascript.js')
copyfile ('./templates/manual.css','./manual_docs/manual.css')

% Update file references
system ('sed -i '''' -e  ''s/octave-forge.css/manual.css/g'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''s/octave-forge.css/manual.css/g'' ./manual_docs/statistics-resampling/function/*.html');

% Create accessible link 
if exist ('./manual.html', 'file')
  delete ('./manual.html')
end
system ('ln -s ./manual_docs/statistics-resampling/index.html ./manual.html');

% Clear uplus
clear myopt
