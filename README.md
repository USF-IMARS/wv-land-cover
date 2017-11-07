# wv2-processing

## filestructure notes
* `./submit.py` is the script which handles command line arguments so we can use scripts from a shell command line.
* `wv2_processing` contains python scripts which do the actual work (the stuff that used to be done by matlab)

## github basics

#### download repo to local machine
`git clone https://github.com/USF-IMARS/wv2-processing`

#### basic git/github workflow
1. `git pull origin master` - this updates your local to match the remote
2. make your file edits
3. `git status` to review the changes you have made
4. (optional) `git diff` to review even more closely
5. `git add my-new-file.py` to add new files to the "staging area"
6. `git commit -a -m "my new commit"` submits a commit with all changes and your commit message "my new commit"
7. `git push origin master` this uploads your commits to github
