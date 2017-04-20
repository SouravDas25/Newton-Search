# Newton-Search
A Enhanced Searching Technique which converge very quickly.
- It uses Newtons inverse interpolation as base algorithms.
- It has a time complexity of O( log <sub>3</sub>( log <sub>3</sub>(n) ) * z ) where z <<< n
- z is the no of column in difference table.
- It uses 3 sub-algorithms to predict the index
  - newtons forward inverse interpolation
  - gaussian middle inverse interpolation
  - newtons backward inverse interpolation
