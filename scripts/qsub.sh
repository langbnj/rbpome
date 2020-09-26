#!/bin/sh

if [[ $@ == "" ]]
then
	echo
	echo "--------------------------------------------------------------------------------";
	echo "> No command given, showing queue instead:";
	echo "--------------------------------------------------------------------------------";
	echo
	/users/gt/blang/scripts/showqueue.pl
	echo
	exit
fi

if [ ! -z $1 ]
then
	# # Absolute path to this script. /home/user/bin/foo.sh
	# DIR=$PWD
	# echo "DIR $DIR"
	# # Absolute path this script is in. /home/user/bin
	# DIR=`dirname $DIR`
	# echo "DIR $DIR"
	# # Strip /users/gt/blang/
	# DIR=`echo $DIR | sed 's/^\/home\/lang\/home\/*//'`
	# echo "DIR $DIR"
	# # Get Name
	# NAME=`echo $@`
	# echo "NAME $NAME"
	# # Combine Path and Name
	# NAME="${DIR}/${NAME}"
	# echo "NAME $NAME"
	# # Strip non-alphanumerics
	# NAME=`echo $NAME | sed 's/[^a-zA-Z0-9]/_/g'`
	# echo "NAME $NAME"
	# Absolute path to this script. /home/user/bin/foo.sh
	ORIGNAME="$PWD/`basename $1`"
	# echo "ORIGNAME $ORIGNAME"
	# Strip /home/lang/home/
	ARGS=$@
	shift
	NAME=$ORIGNAME
	NAME=`echo $NAME | sed 's/^\/lmb\/home\/blang\/*//'`
	NAME=`echo $NAME | sed 's/^\/home\/lang\/home\/*//'`
	NAME=`echo $NAME | sed 's/^\/home\/lang\/*//'`
	NAME=`echo $NAME | sed 's/^\/Users\/lang\/*//'`
	NAME=`echo $NAME | sed 's/^\/g\/scb\/bork\/lang\/*//'`
	NAME=`echo $NAME | sed 's/^\/nfs\/users\/gt\/blang\/*//'`
	NAME=`echo $NAME | sed 's/^\/users\/gt\/blang\/*//'`
	NAME=`echo $NAME $@`
	MYNAME=$NAME
	# echo "NAME $NAME"
	# Strip non-alphanumerics
	NAME=`echo $NAME | sed 's/[^a-zA-Z0-9]/_/g'`
	# echo "NAME $NAME"

	# Shorten to 200ish
	# NAMELEN=${#NAME}
	NAME=${NAME:0:200}
fi

if [[ $ORIGNAME == "" && $ARGS != "" ]]
then
	echo
	echo "--------------------------------------------------------------------------------";
	echo "> Error: '$ARGS' doesn't exist, cancelling";
	echo "--------------------------------------------------------------------------------";
	echo
	exit
fi

if [[ $NAME == "" && $ARGS != "" ]]
then
	echo
	echo "--------------------------------------------------------------------------------";
	echo "> Error: '$ORIGNAME' parsed to '$NAME', cancelling";
	echo "--------------------------------------------------------------------------------";
	echo
	exit
fi

echo "#!/bin/sh" > log-command-$NAME.txt
echo "export START=\`date +%s\`" >> log-command-$NAME.txt
# echo "export TERM=xterm-color" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt
echo "echo \"> $ARGS (\`hostname\`)\"" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "source /users/gt/blang/.bash_profile" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "$ARGS" >> log-command-$NAME.txt
echo "export STOP=\`date +%s\`" >> log-command-$NAME.txt
echo "export ELAPSED=\`expr \$STOP - \$START\`" >> log-command-$NAME.txt
echo "export DAYS=\`expr \$ELAPSED / 86400\`" >> log-command-$NAME.txt
echo "export HOURS=\`expr \$ELAPSED % 86400 / 3600\`" >> log-command-$NAME.txt
echo "export MINS=\`expr \$ELAPSED % 86400 % 3600 / 60\`" >> log-command-$NAME.txt
echo "export SECS=\`expr \$ELAPSED % 86400 % 3600 % 60\`" >> log-command-$NAME.txt
echo "echo" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt
echo "echo \"> Done! (Time elapsed: \$DAYS days \$HOURS h \$MINS min \$SECS sec)\"" >> log-command-$NAME.txt
echo "echo \"================================================================================\"" >> log-command-$NAME.txt

rm -f log-errors-$NAME.txt
rm -f log-output-$NAME.txt
echo 
echo
echo "--------------------------------------------------------------------------------";
echo "> Submitting job '$NAME'";
echo "--------------------------------------------------------------------------------";
echo
qsub -cwd -q long-sl7 -pe smp 1 -l h_rt=720:00:00 -l h_vmem=7G -l virtual_free=7G -N $NAME -e log-errors-$NAME.txt -o log-output-$NAME.txt -S /bin/bash log-command-$NAME.txt
# qsub -cwd -q short-sl65,short-sl7 -pe smp 1 -l h_rt=6:00:00 -l h_vmem=5G -l virtual_free=5G -N $NAME -e log-errors-$NAME.txt -o log-output-$NAME.txt -S /bin/bash log-command-$NAME.txt
# qsub -cwd -q short-sl7 -pe smp 1 -l h_rt=6:00:00 -l h_vmem=7G -l virtual_free=7G -N $NAME -e log-errors-$NAME.txt -o log-output-$NAME.txt -S /bin/bash log-command-$NAME.txt
echo

echo `date +%Y-%m-%d\ %H:%M:%S` "	$ARGS" >> /users/gt/blang/scripts/commandlog.txt
