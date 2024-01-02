% THIS SCRIPT IS ONLY INTENDED FOR USE BY THE PACKAGE MAINTAINER
% Run GNU Octave from a unix-based OS with pandoc installed. 
% Octave needs to have the following packages installed:
%  - generate_html
%  - statistics
%  - statistics-resampling
% (The documentation was generated with GNU Octave installed with MacPorts on MacOS)

% Load required packages
pkg load generate_html
pkg load statistics
pkg load statistics-resampling

% Get options (use Octave-Forge as initial template then customise)
myopt = get_html_options ("octave-forge");

% Customize the options
myopt.include_package_news = false;
myopt.include_alpha = false;
myopt.include_demos = true;
myopt.download_link = "https://gnu-octave.github.io/packages/statistics-resampling/";
myopt.repository_link = "https://github.com/gnu-octave/statistics-resampling/tree/master";
myopt.older_versions_download = "https://github.com/gnu-octave/statistics-resampling/releases";
generate_package_html ('statistics-resampling', 'manual_docs', myopt);

% Delete files that are surplus to requirements
delete ('./manual_docs/footer.js')
delete ('./manual_docs/news.png')
delete ('./manual_docs/manual.png')
delete ('./manual_docs/oct.png')
delete ('./manual_docs/favicon.ico')

% Replace the octave-forge javascript files with custom package templates
delete ('./manual_docs/javascript.js')
delete ('./manual_docs/octave-forge.css')
copyfile ('./templates/javascript.js','./manual_docs/javascript.js')
copyfile ('./templates/manual.css','./manual_docs/manual.css')

% Update file references
system ('sed -i '''' -e  ''s/octave\-forge\.css/manual\.css/g'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''s/octave\-forge\.css/manual\.css/g'' ./manual_docs/statistics-resampling/function/*.html');

% Change page headings
system ('sed -i '''' -e  ''s/statistics\-resampling\<\/h2\>/About this package\<\/h2\>/g'' ./manual_docs/statistics-resampling/index.html');

system ('sed -i '''' -e  ''s/statistics\-resampling\<\/h2\>/Function Reference\<\/h2\>/g'' ./manual_docs/statistics-resampling/overview.html');

% Add left menu to package information page
system ('sed -i '''' -e  ''s/javascript\:fix_top_menu ()\;/javascript\:fix_top_menu ()\; javascript\:show_left_menu ()\;/g'' ./manual_docs/statistics-resampling/index.html');

% Move the icon over
copyfile ('../doc/icon_48x48.png','./manual_docs/icon_48x48.png')

% Change to main package directory
cd ..

% Create accessible link to manual
if exist ('MANUAL.html', 'file')
  delete ('./MANUAL.html')
end
system ('ln -s ./man/manual_docs/statistics-resampling/overview.html MANUAL.html');

% Publish README markdown file as HTML
system ('pandoc README.md -o ./man/manual_docs/statistics-resampling/tmp.html');

% Add header to the readme page
system ('cat ./man/templates/header.html ./man/manual_docs/statistics-resampling/tmp.html > ./man/manual_docs/statistics-resampling/readme.html');
delete ('./man/manual_docs/statistics-resampling/tmp.html')

% Add hyperlink to manual directly in the HTML of the readme file (since the
% link in the markdown file would fail on GitHub page) 
system ('sed -i '''' -e  ''s/Function Reference/\<a href\=\"overview\.html\"\>Function Reference\<\/a\>/g'' ./man/manual_docs/statistics-resampling/readme.html');

% Create accessible link to readme page in manual
if exist ('./README.html', 'file')
  delete ('./README.html')
end
system ('ln -s ./man/manual_docs/statistics-resampling/readme.html README.html');

% Change back to man directory
cd man

% Update meta information in all html files
system ('sed -i '''' -e  ''s/\<meta name\=\"author\" content\=\"\(.*\)\" \/\>/\<meta name\=\"author\" content\=\"Andrew Penn\" \/\>/'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''s/\<meta name\=\"author\" content\=\"\(.*\)\" \/\>/\<meta name\=\"author\" content\=\"Andrew Penn\" \/\>/'' ./manual_docs/statistics-resampling/function/*.html');
system ('sed -i '''' -e  ''s/\<meta name\=\"description\" content\=\"\(.*\)\" \/\>/\<meta name\=\"description\" content\=\"A package for statistical analysis using resampling methods.\" \/\>/'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''s/\<meta name\=\"description\" content\=\"\(.*\)\" \/\>/\<meta name\=\"description\" content\=\"A package for statistical analysis using resampling methods.\" \/\>/'' ./manual_docs/statistics-resampling/function/*.html');
system ('sed -i '''' -e  ''s/\<meta name\=\"keywords\" lang\=\"en\" content\=\"\(.*\)\" \/\>/\<meta name\=\"keywords\" lang\=\"en\" content\=\"GNU Octave Packages, MATLAB Toolbox\" \/\>/'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''s/\<meta name\=\"keywords\" lang\=\"en\" content\=\"\(.*\)\" \/\>/\<meta name\=\"keywords\" lang\=\"en\" content\=\"GNU Octave Packages, MATLAB Toolbox\" \/\>/'' ./manual_docs/statistics-resampling/function/*.html');
system ('sed -i '''' -e  ''/\<link rel\=\"shortcut icon\" href\=\"\.\.\/favicon.ico\" \/\>/d'' ./manual_docs/statistics-resampling/*.html');
system ('sed -i '''' -e  ''/\<link rel\=\"shortcut icon\" href\=\"\.\.\/.\.\/favicon.ico\" \/\>/d'' ./manual_docs/statistics-resampling/function/*.html');

% Clear up
clear myopt
