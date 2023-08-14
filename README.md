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

## Running R Studio Server
- Launch rocker rstudio docker image [source](https://rocker-project.org/)
```bash
# Note: Replace "Data1" with path to whatever folder you want mounted and visible in R
 sudo docker run --rm -ti -v /Data1:/home/rstudio/Data1 -e PASSWORD=jupyter -e ROOT=true -p 8787:8787 rocker/rstudiov2
```
- In browser on local machine: http://54.241.254.117:8787/
- Login with username: rstudio and password entered in command above.

## Running Jypyter Notebook Server
Source: https://docs.aws.amazon.com/dlami/latest/devguide/setup-jupyter-configure-client-windows.html
- Configuring commands:
```bash
# Confugire SSL certificate
openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout mykey.key -out mycert.pem
# Configure jupyter notebook password
jupyter notebook password
```
- Launch server:
```bash
# Launch server (remember to navigate to correct root folder)
jupyter notebook --certfile=~/ssl/mycert.pem --keyfile ~/ssl/mykey.key
```
- Connect to server from local machine:
```bash
# From terminal:
ssh -i "H:\My Drive\AWS\yuri-malina.pem" -N -f -L 8888:localhost:8888 ubuntu@ec2-54-241-254-117.us-west-1.compute.amazonaws.com
# From browser:
https://localhost:8888/
```

### Connecting to SQLite server
source: https://github.com/coleifer/sqlite-web
- Configuring commands
```bash
 pip install sqlite-web
```

- Launch server
```bash
sqlite_web --ad-hoc --password path_to_db 
```
- Connect to server from local machine:
Run the following on local terminal:
```bash
ssh -i "H:\My Drive\AWS\yuri-malina.pem" -N -f -L 8080:localhost:8080 ubuntu@ec2-54-241-254-117.us-west-1.compute.amazonaws.com
```
Note you may need to replace the port number with the appropriate number from the output of the command used to launch the server.

- Connect to server from local machine:
https://localhost:8080/
