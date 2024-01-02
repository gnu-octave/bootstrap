/* Copyright (C) 2008 Soren Hauberg

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
<http://www.gnu.org/licenses/>.
*/

function goto_url (selSelectObject)
{
    if (selSelectObject.options[selSelectObject.selectedIndex].value != "-1") {
        location.href=selSelectObject.options[selSelectObject.selectedIndex].value;
    }
}

function show_left_menu ()
{
  document.getElementById ("left-menu").style.display = "block";
}

function manual_menu ()
{
  // XXX: What should we do here? And do we even need this function?
  write_left_menu ();
}

function write_top_menu (prefix)
{
  // default prefix (maybe some old browsers need this way)
  prefix = (typeof prefix == 'undefined') ? '.' : prefix;

  document.write
  (`
  <div id="top-menu" class="menu"> 
   <table class="menu">
      <tr>
        <td style="width: 90px;" class="menu" rowspan="2">
          <a name="top">
          <img src="${prefix}/oct.png" alt="Octave logo" />
          </a>
        </td>
        <td class="menu" style="padding-top: 0.9em;">
          <big class="menu">Octave Packages</big><small class="menu"> - Extra packages for GNU Octave</small>
        </td>
      </tr>
      <tr>
      </tr>
    </table>
  </div>
   `);
}

function write_docs_left_menu (prefix)
{
  // default prefix (maybe some old browsers need this way)
  prefix = (typeof prefix == 'undefined') ? '.' : prefix;

  document.write
  (`
<div id="left-menu">
  <h3>Navigation</h3>
  <p class="left-menu"><a class="left-menu-link" href="${prefix}/operators.html">Operators and Keywords</a></p>
  <p class="left-menu"><a class="left-menu-link" href="${prefix}/function_list.php">Function List:</a>
  <ul class="left-menu-list">
    <li class="left-menu-list">
      <a  class="left-menu-link" href="${prefix}/octave/overview.html">&#187; Octave core</a>
    </li>
    <li class="left-menu-list">
      <a  class="left-menu-link" href="${prefix}/functions_by_package.php">&#187; by package</a>
    </li>
    <li class="left-menu-list">
      <a  class="left-menu-link" href="${prefix}/functions_by_alpha.php">&#187; alphabetical</a>
    </li>
  </ul>
  </p>
  <p class="left-menu"><a class="left-menu-link" href="${prefix}/doxygen/html">C++ API</a></p>
</div>
   `);
}
