#! /bin/csh

set dir = "."
set N = 204
set x = `ls $dir/A10S-A.job*.log | sort -n -t b -k 2`


set total = `cat $x | grep "EXCHANGE_ACCEPT" | tail -1 | awk '{print $NF/2}'`

set list = `grep -A $N "set cf" $dir/move_pi.conf | tail -$N | tr "{}" " " | awk '{print $1}'`
set klist = `grep -A $N "set cf" $dir/move_pi.conf | tail -$N | tr "{}" " " | awk '{print $2}'`

#echo $list

echo "center_1 \t center2 \t k_1 \t\t k_2 \t accepted_attempts \t total_attempts \t ratio"
set j = 1
while ($j < $N)
set i = `expr $j - 1`
set pre = `echo $list | awk -v i=$i '{print $(i+1)}'`
set nex = `echo $list | awk -v j=$j '{print $(j+1)}'`

set prek = `echo $klist | awk -v i=$i '{print $(i+1)}'`
set nexk = `echo $klist | awk -v j=$j '{print $(j+1)}'`
 
set accept = `cat $x | grep "EXCHANGE_ACCEPT $i $j"  | wc -l`
set r = `echo $accept / $total | bc -l`
echo "$pre \t\t $nex \t\t  $prek \t  $nexk \t\t  $accept \t\t $total \t\t $r"
@ j++
end
