# Goal of these instructions
The aim is to help one build and run 3 particle-in-cell codes (WarpX, Smilei and EPOCH) on a laptop with Linux or on Cineca supercomputers.

# Laptop
## General dependencies
First, install on your laptop the following fundamental software:
```bash
sudo apt-get install git python3-h5py python3-ipython python3-pint python3-sphinx python3-matplotlib python3-dev python3-numpy python3-scipy python3-pip build-essential gcc libhdf5-openmpi-dev 
```
then, follow the related instructions for each case.

To run efficiently in parallel on your machine you need to know your architecture.
For example, you can find out the number of threads per core and cores per socket on your machine with the command: `lscpu`

# Cineca
You must have a Cineca account and two-factor authentication (2FA) enabled (https://wiki.u-gov.it/confluence/display/SCAIUS/How+to+connect+via+2FA).

To access Galileo100:
```bash
ssh <username>@login.g100.cineca.it
```
similar for Marconi100:
```bash
ssh <username>@login.m100.cineca.it
```
and Leonardo:
```bash
ssh <username>@login.leonardo.cineca.it
```
then follow the related instructions. 

# Git 
There are many tutorials online, please consider this just a summary/reminder of the main git commands.

If you want to contribute to the `PIC_tools` repo, first fork it in your GitHub account, then clone it: 
```bash
git clone git@github.com:<your_git_username>/PIC_tools.git
```
to incorporate changes from a remote repository into the current branch: 
```bash
git pull
```
to prepare the content staged for the next commit:
```bash
git add <files_to_add>
```
create a new commit, i.e. record the latest changes of the source code to your repository:
```bash
git commit -a -m "useful description"
```
update remote repository:
```bash
git push
```

If you want to merge your changes in the main repository, go to GitHub and make a pull request:
If the pull request is accepted and the code is merged, then go back to your repo and synch your fork, then do git pull from the terminal to update your remote repository. 
