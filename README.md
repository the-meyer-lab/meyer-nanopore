# Meyer Lab Nanopore Project
## Purpose  
The purpose of this repository is to create a data pipeline for the processing and analysis of nanopore sequence data.

## Setup
### 1. Configure AWS accounts
*(skip if you've already configured your AWS account and key pair)*
#### a. Set up new admin account (if applicable)
- Organizations > add an aws account
- Follow instructions

#### b. Set up a new IAM user for individual users
- IAM > users > add users
- Follow instructions

#### c. Create key pair for users who will be accessing instances
*(This is only required once for each user, and must be kept private.*)
- EC2 > Key pairs >
	- Enter name of user key pair will be granted to
	- Leave defaults: RSA / .pem
	- Download file, and securely provide to user who will be using it. Delete it once transfer is complete.
	
### 2. Configure AWS services
*(These steps can be skipped if your EC2 instances, EBV and S3 buckets are already configured.)*
#### a. Configure a new IAM role that will be granted to EC2 instances allowing for S3 access:
*(This is only required once when setting up AWS account for first time.)*
- IAM > Create role > AWS Service
	- Select use case, in this case allowing EC2 to call AWS services on our behalf.
- Permissions:
	- Select 'AmazonS3FullAccess'
- Name role 'EC2toS3' 
- Create role

#### b. Set up Elastic Block Volume (EBV, hard drive that will be attached to our compute instance)
- EC2 > Elastic Block Store > Volumes > Create Volume
- Select volume type (General Purpose SSD is suitable for most applications)
- Select desired size
- Select region (needs to match region where you will be using instance.) For Meyer Lab applications, select us-west-1b

#### c. Set up compute instance
- If desired instance is already configured, select from list and launch. If notm configure new instance with following steps:
	- EC2 > Instances > Launch Instance > 
		- Select operating system: recommend Deep Learning AMI (Ubuntu 18.XX) Version XX.XX
		- Select instance type: consult on-demand instance pricing [here](https://aws.amazon.com/ec2/pricing/on-demand/)
			- For compute intensive workloads (e.g. base calling), recommend: g4dn.metal @ ~$9/hr
			>Note: When configuring a g4dn.metal instance for the first time, AWS default commercial accounts have a quota of 0 for these types of instances which will cause an error when attempting to launch the instance. In order to resolve this, visit Service Quotas > EC2 > request quota increase for "All G and VT Spot Instance Requests" up to '96' (which is the amount required by a g4dn.metal instance). 
			- For less intensive general purpose workloads (e.g. alignment), recommend:  t3a.2xlarge @ 	$0.39/hr
		- Configure instance details: 
			- Pick subnet that matches the region of the EBV. For Meyer Lab applications select us-west-1b
			- Assign IAM role "EC2toS3"
		- Review & Launch > Launch
			- Assign key pair of users who will be accessing instance. Create new key pair if required by following [[Nanopore Sequencing AWS Setup Manual#c Create key pair for users who will be accessing instances]]

#### d. Attach EBV to instance
- EBS > Volumes
	- Select volume created in [[Nanopore Sequencing AWS Setup Manual#b Set up Elastic Block Volume EBV hard drive that will be attached to our compute instance]]
	- Actions > Attach Volume >
		- Select instance

#### e. Set up S3 bucket
- S3 > Create bucket
	- Add bucket name. Recommend a new bucket for each project.
	- Select region, ensure matches previously selected regions of instances & EBV.
	- Select whether versioning is desired. Disabled by default is fine.
		-  See standard S3 pricing [here](https://aws.amazon.com/s3/pricing/) to help you decide whether to use versioning. Versioning greatly increases storage utilization, see explanation [here](https://aws.amazon.com/s3/faqs/).
	- Create bucket.

### 3. Uploading your data to s3
- From the AWS console, navigate to your S3 bucket created in [[Nanopore Sequencing AWS Setup Manual#e Set up S3 bucket]]
- Click "Upload" and "Upload Folder"
- Select folder containing your raw fast5 folders from nanopore sequencing.

### 4. Configuring your instance
#### a. Install AWS CLI
- Follow instructions on installing AWS CLI (Amazon Commandline Interface) [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

#### b. Connect to instance from your command line
- EC2 > Instances > select desired instance that has been launched and click "Connect"
	- Will display a sample command to connect to your instance e.g.:
	```bash
	ssh -i "G:\My Drive\Meyer Lab\AWS\yuri-malina.pem" ubuntu@ec2-54-219-137-86.us-west-1.compute.amazonaws.com
	```
	- "G:\My Drive\Meyer Lab\AWS\yuri-malina.pem" should be replaced by the path to your .pem file downloaded in step [[Nanopore Sequencing AWS Setup Manual#c Create key pair for users who will be accessing instances]].
		- "c2-54-219-137-86" should be replaced with the address listed on your instance's connect page.
	- If you get an error with this command, an alternative is:
	```bash
	ssh -i "G:\My Drive\Meyer Lab\AWS\yuri-malina.pem" ec2-c2-54-219-137-86.us-west-1.compute.amazonaws.com -l ubuntu
	```

#### c. Configure your EVB drive on your AWS instance
Excute the following commands:

```bash
###### Configure EBV drive on instance  
# Format drives to the correct file system type, and mount them  
# Note: only do this when mounting drive for first time. This will clear data from an existing drive.
# For drives < 2TB in size:
printf "n\np\n\n\n\nw" | sudo fdisk /dev/nvme1n1  
# For drives > 2TB in size:
sudo mkfs -t xfs /dev/xvdf
# fdisk options: n = new partition; p = primary; accept 3 defaults; w = write  
sudo partprobe /dev/nvme1n1  
printf "y" | sudo mkfs.ext4 /dev/nvme1n1

# For copying large amount of files between drives: rsync -ah --info=progress2 [source] [destination]  
  
# Make folder (if it doesn't already exist) and mount drive  
sudo mkdir -p /Data1  
# Determine your drive name (usually of the format nvmeXnX)
# to list available drives: sudo lsblk
sudo mount /dev/nvme1n1 Data1  
# Change owner of drive  
sudo chown -R ubuntu Data1  
######  
```

#### d. Clone this git repository onto the instance
Execute the following commands if installing on ubuntu ([source](https://github.com/cli/cli/blob/trunk/docs/install_linux.md)):
- Install github command line
```bash
mkdir -p /Data1/git  
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null
sudo apt update
sudo apt install gh
```
For other OS see gh CLI installation instructions [here](https://github.com/cli/cli#installation).
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

#### e. Install libraries and assets required for nanopore sequencing 
Execute python script from this git repo, by navigating to your script 

```bash
sh ../git/meyer-nanopore/scripts/instance_initial_config.sh
```

### 5. Quality of life (optional tools)
When working with a remote server, if you are not entirely comfortable with the command line interface, navigating around, manipulating files, and creating / editing text files with vim, nano or other CLI text editor, it is highly recommend to install some tools to improve your workflow:
- SFTP / FTP Client: *Allows you to view your local and remote file structure in a GUI, allowing for easy drag and drop file transfers.*
	- Windows:  [WinSCP](https://winscp.net/eng/index.php) is a good choice but there are others. 
	- Mac: [Cyberduck](https://cyberduck.io/)
- Integrated Development Environment: *Gives you a fully featured code editor, console window and file structure all in one interface*
	- Windows or Mac: I recommend [JetBrains Gateway](https://www.jetbrains.com/remote-development/gateway/) (free for .edu email addresses) but again there are others.

## Usage
*...To follow...*


{
   "apiKey": "",
   "project_id": "igv",

   "auth_provider": "Amazon",
   "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
   "aws_region": "ap-southeast-2",
   "scope": "email%20openid%20profile",
   "redirect_uris": [
     "http://localhost:60151/oauthCallback"
   "client_id": "7imktrjfl5p1to6kbiceo115pk",
   "client_secret": "1ign75opko356j563gdsh57k8rkfkgkdv0l9b115ir9lkm005omk",
   "authorization_endpoint": "https://meyer-lab-igv.auth.us-west-1.amazoncognito.com/login",
   "token_endpoint": "https://meyer-lab-igv.auth.us-west-1.amazoncognito.com/token",
   "aws_cognito_fed_pool_id": "us-west-1:71c57749-c5b5-4b58-8034-b70195275c09",
   "aws_cognito_pool_id": "us-west-1_Bv21MeV6i",
   "aws_cognito_role_arn": "arn:aws:iam::995302141357:role/Cognito_igvidpoolAuth_Role"
}
