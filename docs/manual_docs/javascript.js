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
        <td class="menu" style="padding-top: 0em;">
          <big class="menu">Manual</big>
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
  <p class="left-menu"><a class="left-menu-link" href="${prefix}/statistics-resampling/index.html">Package Information</a></p>
  <p class="left-menu"><a class="left-menu-link" href="${prefix}/statistics-resampling/overview.html">Function Reference</a></p>
</div>
   `);
}
