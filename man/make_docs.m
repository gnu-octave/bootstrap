% THIS SCRIPT IS ONLY INTENDED FOR USE BY THE PACKAGE MAINTAINER
% Run GNU Octave from a Bash terminal in a unix-based OS with pandoc installed. 
% Octave needs to have the following packages installed:
%  - generate_html
%  - statistics
%  - statistics-resampling
% (The documentation was generated with GNU Octave installed with MacPorts on MacOS)

% Remove existing pages and files
system ('rm -Rf ../docs/*');
system ('rm -Rf ../tmp');

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
myopt.repository_link = "https://github.com/gnu-octave/statistics-resampling/";
myopt.older_versions_download = "https://github.com/gnu-octave/statistics-resampling/releases";
myopt.index_body_command = myopt.body_command;
myopt.overview_filename = 'function_reference.html';
generate_package_html ('statistics-resampling', 'tmp', myopt);

% Run (bash) shell script
system ('./make_docs_helper.sh');

% Clear up
clear myopt
