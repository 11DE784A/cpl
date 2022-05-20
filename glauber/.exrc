if &cp | set nocp | endif
let s:cpo_save=&cpo
set cpo&vim
inoremap <silent> <C-Tab> =UltiSnips#ListSnippets()
imap <C-G>S <Plug>ISurround
imap <C-G>s <Plug>Isurround
imap <C-S> <Plug>Isurround
snoremap <silent>  "_c
nnoremap  
xnoremap <silent> 	 :call UltiSnips#SaveLastVisualSelection()gvs
snoremap <silent> 	 :call UltiSnips#ExpandSnippetOrJump()
nnoremap <NL> <NL>
nnoremap  
nnoremap  
snoremap  "_c
xmap S <Plug>VSurround
nmap cS <Plug>CSurround
nmap cs <Plug>Csurround
nmap ds <Plug>Dsurround
xmap gx <Plug>NetrwBrowseXVis
nmap gx <Plug>NetrwBrowseX
xmap gS <Plug>VgSurround
nmap ySS <Plug>YSsurround
nmap ySs <Plug>YSsurround
nmap yss <Plug>Yssurround
nmap yS <Plug>YSurround
nmap ys <Plug>Ysurround
xnoremap <silent> <Plug>NetrwBrowseXVis :call netrw#BrowseXVis()
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#BrowseX(netrw#GX(),netrw#CheckIfRemote(netrw#GX()))
snoremap <C-R> "_c
snoremap <silent> <C-H> "_c
snoremap <silent> <Del> "_c
snoremap <silent> <BS> "_c
snoremap <silent> <C-Tab> :call UltiSnips#ListSnippets()
nnoremap <silent> <Plug>SurroundRepeat .
nnoremap <C-H> 
nnoremap <C-L> 
nnoremap <C-K> 
nnoremap <C-J> <NL>
imap S <Plug>ISurround
imap s <Plug>Isurround
inoremap <silent> 	 =UltiSnips#ExpandSnippetOrJump()
imap  <Plug>Isurround
inoremap jk 
let &cpo=s:cpo_save
unlet s:cpo_save
set background=dark
set backspace=indent,eol,start
set backupdir=~/.cache/vim/backup//
set confirm
set directory=~/.cache/vim/swap//
set fileencodings=ucs-bom,utf-8,default,latin1
set helplang=en
set hidden
set ignorecase
set laststatus=2
set ruler
set runtimepath=~/.vim,~/.vim/plugged/gruvbox,~/.vim/plugged/vim-surround,~/.vim/plugged/vim-gutentags,~/.vim/plugged/julia-vim,~/.vim/plugged/vim-cpp-modern,~/.vim/plugged/ultisnips,/usr/share/vim/vimfiles,/usr/share/vim/vim82,/usr/share/vim/vimfiles/after,~/.vim/plugged/vim-cpp-modern/after,~/.vim/plugged/ultisnips/after,~/.vim/after
set shiftwidth=4
set showcmd
set smartcase
set smarttab
set spelllang=en_gb
set splitbelow
set splitright
set statusline=%1*\ buf:%n\ %2*%{b:gitstatus}%*\ %F%m%r%h%w%=%y\ %2*\ %{&ff}\ %1*\ ln:%l/%L\ %p%%\ 
set suffixes=.bak,~,.o,.info,.swp,.aux,.bbl,.blg,.brf,.cb,.dvi,.idx,.ilg,.ind,.inx,.jpg,.log,.out,.png,.toc
set tabstop=4
set termguicolors
set undodir=~/.cache/vim/undo//
" vim: set ft=vim :
