# TMSMultiLab
A place to share resources, data and code for running a transcranial magnetic stimulation (TMS) research laboratory...

## GitHub tutorial
If you've never used git or github before, you'll need some tutorials:

https://docs.github.com/en/get-started/start-your-journey/hello-world

If you're happy to get started, then there are three main ways you can use the repository.

## 1. Download the files you need
If you just want a few files, navigate to them and find the download button. Enjoy!

## 2. Set up a read-only copy of the TMSMultiLab repository
You can use git commands (or the github website) to clone the whole repository. This will set up all the TMSMultiLab files and directories on your local machine.

Later on, when you want to update your local version to the latest versions, you can run the code <code>git pull</code> from the TMSMultiLab directory on your local machine and it should do everything you need. Here's a quick guide to cloning the repository and the wiki:

### Cloning the repository to your local machine

1. install git on your local machine
2. open a terminal
3. navigate to the directory where you want the repository to be cloned (eg: MyFiles)
4. type <code>git clone  https://github.com/TMSMultiLab/TMSMultiLab</code>
5. it is done - you now have all the code in <code>MyFiles/TMSMultiLab</code>

If you edit these files on your local machine, you will then need to <code>stage</code>, <code>commit</code>, and <code>push</code> them up to the main repository - this needs to be done as per below = 3. Developing the TMSMultiLab code & wiki

To ignore any local changes when trying to pull updated files from the repository, you can use <code>git checkout Directory/File.ext</code> to ignore local changes on each file, and overwrite this local file with the updated repository file.

If the main repository files change, you will need to <code>fetch</code> or <code>pull</code> the updates before using them.

### Cloning the wiki to your local machine
see 1-3 above, then:

4. type <code>git clone https://github.com/TMSMultiLab/TMSMultiLab.wiki.git</code>
5. it is done - you now have all the wiki in <code>MyFiles/TMSMultiLab.wiki</code>

If the wiki files change, you will need to <code>fetch</code> or <code>pull</code> the updates before working on them.

## 3. Developing the TMSMultiLab code & wiki
If you're planning to make major changes to the code, or to contribute new material, you will most likely want to create your own <code>fork</code> of the TMSMultiLab repository. This allows you to work on the files without affecting anyone else's version of the repository. Once your changes have been completed, tested, and reviewed with other TMSMultiLab members, they can then be merged across into the main repository.

1. Log in to github and browse to the [TMSMultiLab repository](https://github.com/TMSMultiLab/TMSMultiLab)
2. Click on the <code>Fork</code> menu button, choose the account to which you want the repository to be forked, and include the <code>main</code> branch
3. You can now edit your own complete version of the TMSMultiLab code and wiki. While you _could_ do this on the github website, most likely you'll do this in a cloned version of your own (forked) repository on your local machine
4. Come back here when you're ready to share and merge your changes...
