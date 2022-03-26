# Meyer Lab Nanopore Sequencing
## Purpose  
The purpose of this repository is to create a data pipeline for the processing and analysis of nanopore sequencing datasets.

## Requirements
### 1. AMI
These scripts must be executed on an AWS instance launched with the following AMI:
ami-07ac92c40cd437ed0

See detailed in packages included in this AMI in this github repository:
https://github.com/the-meyer-lab/aws-cloud-computing-setup

If instance was not launched form this AMI, you must execute all commands in the "nanopore-analysis-ami.sh" in the above repository.

### 2. Clone this git repository onto the instance
- Configure git with your git user credentials. If you don't have git credentials, follow instructions [here](https://docs.github.com/en/get-started/signing-up-for-github/signing-up-for-a-new-github-account) to create them for the first time. Replace "frnkmxwll" and "ymalina@gmail.com" with your user name and email used to create your github account.
```bash
git config --global user.name "frnkmxwll"
git config --global user.email ymalina@gmail.com
```
- Authenticate using Github Commandline
```bash
gh auth login
# Log into: **GitHub.com**
# Preferred protocol is: HTTPS
# Authenticate Git with GitHub credentials: Yes
# How would you like to authenticate: Paste an authentication token
```

To obtain an authentication token follow github instructions [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
	Select the following scopes:
	![[Pasted image 20220121191329.png]]

Clone this git repository with the following command:

NOTE: only necessary if using new EBV volume and git repository is not already cloned. 
If repository already exists, execute "git pull" to update local branch.
```bash
git clone https://github.com/frnkmxwll/meyer-nanopore ../git/

git pull
```

## Basecalling and analysis
Execute commands in "basecalling.sh" file.

More detailed instructions to follow...