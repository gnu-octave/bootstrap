% This script is only intended for use by the package maintainer
pkg load generate_html
pkg load parallel
pkg load statistics
pkg load statistics-resampling
myopt = get_html_options ("octave-forge");
myopt.include_package_news = false;
myopt.download_link = "https://gnu-octave.github.io/packages/statistics-resampling/";
myopt.repository_link = "https://github.com/gnu-octave/statistics-resampling/tree/master";
myopt.older_versions_download = "https://github.com/gnu-octave/statistics-resampling/releases";
generate_package_html ('statistics-resampling', 'package_documentation', myopt);
