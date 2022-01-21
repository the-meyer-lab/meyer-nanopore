# Meyer Lab Nanopore Project  
## Purpose  
The purpose of this repository is to create a data pipeline for the processing and analysis of nanopore sequence data.

## Usage

## Setup
## 1. Configure AWS accounts
*(skip if you've already configured your AWS account and key pair)*
### a. Set up new admin account (if applicable)
- Organizations > add an aws account
- Follow instructions

### b. Set up a new IAM user for individual users
- IAM > users > add users
- Follow instructions

### c. Create key pair for users who will be accessing instances
*(This is only required once for each user, and must be kept private.*)
- EC2 > Key pairs >
	- Enter name of user key pair will be granted to
	- Leave defaults: RSA / .pem
	- Download file, and securely provide to user who will be using it. Delete it once transfer is complete.
	
## 2. Configure AWS services
*(These steps can be skipped if your EC2 instances, EBV and S3 buckets are already configured.)*
### a. Configure a new IAM role that will be granted to EC2 instances allowing for S3 access:
*(This is only required once when setting up AWS account for first time.)*
- IAM > Create role > AWS Service
	- Select use case, in this case allowing EC2 to call AWS services on our behalf.
- Permissions:
	- Select 'AmazonS3FullAccess'
- Name role 'EC2toS3' 
- Create role

### b. Set up Elastic Block Volume (EBV, hard drive that will be attached to our compute instance)
- EC2 > Elastic Block Store > Volumes > Create Volume
- Select volume type (General Purpose SSD is suitable for most applications)
- Select desired size
- Select region (needs to match region where you will be using instance.) For Meyer Lab applications, select us-west-1b

### c. Set up compute instance
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

### d. Attach EBV to instance
- EBS > Volumes
	- Select volume created in [[Nanopore Sequencing AWS Setup Manual#b Set up Elastic Block Volume EBV hard drive that will be attached to our compute instance]]
	- Actions > Attach Volume >
		- Select instance

### e. Set up S3 bucket
- S3 > Create bucket
	- Add bucket name. Recommend a new bucket for each project.
	- Select region, ensure matches previously selected regions of instances & EBV.
	- Select whether versioning is desired. Disabled by default is fine.
		-  See standard S3 pricing [here](https://aws.amazon.com/s3/pricing/) to help you decide whether to use versioning. Versioning greatly increases storage utilization, see explanation [here](https://aws.amazon.com/s3/faqs/).
	- Create bucket.

## 3. Configuring your instance
### a. Install AWS CLI
- Follow instructions on installing AWS CLI (Amazon Commandline Interface) [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

### b. Connect to instance from your command line
- EC2 > Instances > select desired instance and click "Connect"
	- Will display a sample command to connect to your instance e.g.:
		>ssh -i "yuri-malina.pem" ubuntu@ec2-54-219-34-238.us-west-1.compute.amazonaws.com
	- "yuri-malina.pem" should be replaced by the path to your .pem file downloaded in step [[Nanopore Sequencing AWS Setup Manual#c Create key pair for users who will be accessing instances]].
	- If you get an error with this command, an alternative is:
		> ssh -i "yuri-malina.pem" ec2-54-219-34-238.us-west-1.compute.amazonaws.com -l ubuntu
