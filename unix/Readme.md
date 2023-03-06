

[awk](https://www.geeksforgeeks.org/awk-command-unixlinux-examples/)



```
$ awk '{print}' employee.txt

ajay manager account 45000
sunil clerk account 25000
varun manager sales 50000
amit manager account 47000
tarun peon sales 15000
deepak clerk sales 23000
sunil peon sales 13000
satvik director purchase 80000 
```

```
$ awk '/manager/ {print}' employee.txt 

ajay manager account 45000
varun manager sales 50000
amit manager account 47000
```

```
$ awk '{print $1, "-" , $4}' employee.txt

ajay - 45000
sunil - 25000
varun - 50000
amit - 47000
tarun - 15000
deepak - 23000
sunil - 13000
satvik - 80000
```


| NR: line number; NF: last column

```
$ awk '/account/ {print NR,$0}' employee.txt

1 ajay manager account 45000
2 sunil clerk account 25000
4 amit manager account 47000
```

```
$ awk '/manager/ {print $1, $NF}' employee.txt

ajay 45000
varun 50000
amit 47000
```






## The number of columns

``` head -1 ./data/merge0422_ped-map.ped | tr ' ' '\n' | wc -l ```




## Manipulate column

``` sed 's/_[A-Z].* / /g' file ```


