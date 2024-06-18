#!/usr/bin/env bash

BIN_DIR=bin
CONFIG_DIR=config
SCRIPT_DIR=scripts

programs=( "$BIN_DIR"/* )
configs=( "$CONFIG_DIR"/* )
prog=""
conf=""

# Searches specified config file for the line containing the key string
# then cuts the line after the =
# then removes any leading whitespace and any double-quotation marks.
get_field () {
    local config="$1"
    local key="$2"
    val=$(cat "$config" | grep "$key" | cut -f2 --delimiter='=' | awk '{$1=$1};1' | tr -d '"')
    echo "$val"
}

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
            job_name=$(get_field "$conf" name)
            cpus=$(get_field "$conf" num-cpus)
            max_time=$(get_field "$conf" max-compute-time)
            mem=$((256 / $cpus))

			echo "Scheduling $prog $conf"
			sbatch -J $job_name --ntasks-per-node="$cpus" \
                --time="$max_time" --mem-per-cpu="$mem"G \
                -o slurm-out/"$job_name".out -e slurm-out/"$job_name".err \
                $SCRIPT_DIR/run_job.sh $prog $conf $cpus
			break
			;;
	esac
done

echo "Scheduling done."

