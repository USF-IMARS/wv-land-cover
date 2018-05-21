# wv2-processing
Processing scripts for decision-tree land use classification on WorldView-2 images.
 
The submit_py.sh file is what I use in Circe to call the pgc_ortho.py script, which has a number of sub-scripts called. 
The submit_py.sh also contains the Matlab script call, so you'll want to comment out those lines before testing it.

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
