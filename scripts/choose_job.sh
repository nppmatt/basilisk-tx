#!/usr/bin/env bash

BIN_DIR=bin
CONFIG_DIR=config
SCRIPT_DIR=scripts

programs=( "$BIN_DIR"/* )
configs=( "$CONFIG_DIR"/* )
prog=""
conf=""
select option in "${programs[@]}" "cancel"; do
	case $option in
		"cancel")
            echo "Cancelling job selection."
			exit 0
			;;
		*)
            prog="$option"
			echo "Selected $option"
			break
			;;
	esac
done

select option in "${configs[@]}" "cancel"; do
	case $option in
		"cancel")
            echo "Cancelling job selection."
			exit 0
			;;
		*)
            conf="$option"
			echo "Scheduling $prog $conf"
			sbatch $SCRIPT_DIR/run_job.sh $prog $conf
			break
			;;
	esac
done

echo "Scheduling done."

