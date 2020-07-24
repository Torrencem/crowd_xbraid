let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Documents/Projects/crowd_xbraid
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +212 src/crowd_horesh.cpp
badd +1 MATLAB/main.m
badd +1 term://.//69526:/bin/bash
argglobal
%argdel
$argadd src/crowd_horesh.cpp
edit src/crowd_horesh.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 96 + 119) / 238)
exe '2resize ' . ((&lines * 30 + 28) / 57)
exe 'vert 2resize ' . ((&columns * 141 + 119) / 238)
exe '3resize ' . ((&lines * 23 + 28) / 57)
exe 'vert 3resize ' . ((&columns * 141 + 119) / 238)
argglobal
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let s:l = 387 - ((28 * winheight(0) + 27) / 54)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
387
normal! 08|
wincmd w
argglobal
if bufexists("MATLAB/main.m") | buffer MATLAB/main.m | else | edit MATLAB/main.m | endif
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let s:l = 22 - ((21 * winheight(0) + 15) / 30)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
22
normal! 026|
wincmd w
argglobal
if bufexists("term://.//69526:/bin/bash") | buffer term://.//69526:/bin/bash | else | edit term://.//69526:/bin/bash | endif
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 10023 - ((22 * winheight(0) + 11) / 23)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
10023
normal! 048|
wincmd w
exe 'vert 1resize ' . ((&columns * 96 + 119) / 238)
exe '2resize ' . ((&lines * 30 + 28) / 57)
exe 'vert 2resize ' . ((&columns * 141 + 119) / 238)
exe '3resize ' . ((&lines * 23 + 28) / 57)
exe 'vert 3resize ' . ((&columns * 141 + 119) / 238)
tabnext 1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToOFc
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
