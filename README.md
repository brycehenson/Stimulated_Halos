# Stimulated_Halos
**Bryce M. Henson**  
**Status:** This Code is **in dev and not ready to be used by others**. Unit Testing is **not** implemented.

## To Do
- [ ] Figure out how to share common non packaged utils
  - dont want to have a seperate repo for each
  - could make a basic analysis repo and dump it all in there
- [ ] rebuild functionality with the current tooling
- [ ] build a correlation repo (look at davids code ) with unit tests  


## Install
``` 
git clone --recursive https://github.com/brycehenson/Stimulated_Halos.git
```
then to update 
```
git submodule update --recursive
git submodule foreach --recursive git pull origin master
```



## Contributions  
This project would not have been possible without the many open source tools that it is based on. In no particular order: 

* ***James Conder*** [gaussfilt](https://au.mathworks.com/matlabcentral/fileexchange/43182-gaussfilt-t-z-sigma)
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Jan*** [FileTime](https://au.mathworks.com/matlabcentral/fileexchange/24671-filetime)
* ***Benjamin Kraus*** [nanconv](https://au.mathworks.com/matlabcentral/fileexchange/41961-nanconv)
* ***M. A. Hopcroft**** [allan](https://au.mathworks.com/matlabcentral/fileexchange/13246-allan)
* ***Daniel Eaton***  [sfigure](https://au.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure)
* ***Denis Gilbert***  [M-file Header Template](https://au.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template)
* ***DrosteEffect***  [CIECAM02](https://github.com/DrosteEffect/CIECAM02)
