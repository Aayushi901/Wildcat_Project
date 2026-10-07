#Both PSMC inputs passed validation. We can now run PSMC separately for the two cats, in parallel.
#All 18 autosomes were retained. Missing bins are 42.59% for the flat-headed cat and 27.65% for the rusty-spotted cat. These are valid input files, but validation does not establish that the demographic estimates will be unbiased.
#1. Start both PSMC runs
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

set -o pipefail

printf '%s\n' flat_headed rusty_spotted \
  | xargs -P 2 -n 1 bash -c '
set -euo pipefail

SAMPLE="$1"

case "$SAMPLE" in
    flat_headed) MAXDP=36 ;;
    rusty_spotted) MAXDP=27 ;;
    *) exit 1 ;;
esac

DIR="$SAMPLE/psmc_d10_D${MAXDP}"
INPUT="$DIR/input.psmcfa"
OUT="$DIR/main.psmc"

test -s "$INPUT"

if [ -e "$OUT" ]; then
    echo "STOP: $OUT already exists; preserving it." >&2
    exit 1
fi

echo "START PSMC: $SAMPLE"
date

psmc \
    -N30 \
    -t15 \
    -r5 \
    -p "4+25*2+4+6" \
    -o "$OUT.partial" \
    "$INPUT" \
    > "$DIR/main.stdout.log" \
    2> "$DIR/main.stderr.log"

awk "
    /^RD[[:space:]]/ { round=\$2 }
    /^\/\/$/ { completed=round }
    END {
        if (completed!=30) {
            print \"ERROR: final completed round is not 30\" > \"/dev/stderr\"
            exit 1
        }
    }
" "$OUT.partial"

mv "$OUT.partial" "$OUT"

echo "DONE PSMC: $SAMPLE"
date
' _ \
  && echo "PASS: both PSMC runs completed 30 iterations"

  #2. Optional: check progress from a second terminal
  cd /home/aayushi/data1/psmc_linux

tail -n 10 flat_headed/psmc_d10_D36/main.stderr.log

tail -n 10 rusty_spotted/psmc_d10_D27/main.stderr.log

#3. After the final PASS, print the fitting summaries
cd /home/aayushi/data1/psmc_linux

for SAMPLE in flat_headed rusty_spotted
do
    case "$SAMPLE" in
        flat_headed) DIR="$SAMPLE/psmc_d10_D36" ;;
        rusty_spotted) DIR="$SAMPLE/psmc_d10_D27" ;;
    esac

    echo "=== $SAMPLE ==="

    awk '
    $1=="RD" { round=$2 }
    $1=="LK" || $1=="LL" {
        print round,$2
    }
    ' "$DIR/main.psmc" > "$DIR/likelihood.tsv"

    echo "Last 10 iteration likelihoods:"
    tail -n 10 "$DIR/likelihood.tsv"

    echo "Final iteration and fitted parameters:"
    awk '
    $1=="RD" { rd=$0 }
    $1=="TR" { tr=$0 }
    END {
        print rd
        print tr
    }
    ' "$DIR/main.psmc"
done

#First, run this quick diagnostic on the existing results. It checks:
#- Changes in fitted parameters between iterations 25 and 30.
#- Expected recombination support for each of the 28 fitted time groups.
cd /home/aayushi/data1/psmc_linux

for SAMPLE in flat_headed rusty_spotted
do
    case "$SAMPLE" in
        flat_headed) DIR="$SAMPLE/psmc_d10_D36" ;;
        rusty_spotted) DIR="$SAMPLE/psmc_d10_D27" ;;
    esac

    echo "=== $SAMPLE ==="

    awk '
    $1=="RD" { round=$2 }

    $1=="TR" && (round==25 || round==30) {
        theta[round]=$2
        rho[round]=$3
    }

    $1=="RS" && round==25 {
        previous_lambda[$2]=$4
    }

    $1=="RS" && round==30 {
        i=$2
        current_lambda[i]=$4
        rows++

        if (i<4) group=1
        else if (i<54) group=2+int((i-4)/2)
        else if (i<58) group=27
        else group=28

        support[group]+=$5
    }

    END {
        if (rows!=64 || theta[25]<=0 || theta[30]<=0 ||
            rho[25]<=0 || rho[30]<=0) {
            print "ERROR: missing or unexpected fitting records"
            exit 1
        }

        printf "Theta: iteration 25 = %.6f; iteration 30 = %.6f\n",
               theta[25],theta[30]
        printf "Rho:   iteration 25 = %.6f; iteration 30 = %.6f\n",
               rho[25],rho[30]

        for (i=0; i<64; i++) {
            if (previous_lambda[i]<=0 || current_lambda[i]<=0) {
                print "ERROR: invalid population-size parameter"
                exit 1
            }

            old=theta[25]*previous_lambda[i]
            new=theta[30]*current_lambda[i]
            change=100*(new-old)/old
            if (change<0) change=-change

            if (change>maximum) {
                maximum=change
                worst=i
            }
        }

        printf "Largest relative change in theta*lambda, iterations 25-30: %.2f%%\n",
               maximum
        print "Atomic interval with largest change:",worst

        print "Group\tExpected_recombinations\tAssessment"
        for (g=1; g<=28; g++) {
            status=(support[g]<10 ? "LOW_SUPPORT" : ">=10")
            printf "%d\t%.3f\t%s\n",g,support[g],status
            if (support[g]<10) weak++
        }

        print "Groups below 10 expected recombinations:",weak+0
    }
    ' "$DIR/main.psmc" > "$DIR/convergence_support_round30.txt"

    cat "$DIR/convergence_support_round30.txt"
done

<!--
#results
#(=== flat_headed ===
#Theta: iteration 25 = 0.021289; iteration 30 = 0.020987
#Rho:   iteration 25 = 0.005240; iteration 30 = 0.005279
#Largest relative change in theta*lambda, iterations 25-30: 14.49%
#Atomic interval with largest change: 4
```
Group   Expected_recombinations Assessment
1       23633.382       >=10
2       1764.454        >=10
3       3839.636        >=10
4       4974.857        >=10
5       4678.262        >=10
6       3996.722        >=10
7       3542.419        >=10
8       3657.638        >=10
9       4075.671        >=10
10      4592.785        >=10
11      5141.741        >=10
12      5687.952        >=10
13      6150.512        >=10
14      6493.802        >=10
15      6734.745        >=10
16      6851.958        >=10
17      6776.050        >=10
18      6401.971        >=10
19      5655.229        >=10
20      4624.740        >=10
21      3535.557        >=10
22      2584.882        >=10
23      1855.491        >=10
24      1334.871        >=10
25      970.555 >=10
26      705.119 >=10
27      820.392 >=10
28      424.073 >=10
Groups below 10 expected recombinations: 0
=== rusty_spotted ===
Theta: iteration 25 = 0.037690; iteration 30 = 0.035357
Rho:   iteration 25 = 0.010878; iteration 30 = 0.010689
Largest relative change in theta*lambda, iterations 25-30: 25.20%
Atomic interval with largest change: 4
Group   Expected_recombinations Assessment
1       69549.626       >=10
2       3010.237        >=10
3       9986.415        >=10
4       9838.413        >=10
5       7328.193        >=10
6       5558.381        >=10
7       4772.844        >=10
8       4881.917        >=10
9       5591.730        >=10
10      6807.513        >=10
11      8380.678        >=10
12      10095.963       >=10
13      12152.732       >=10
14      14866.766       >=10
15      18273.180       >=10
16      22104.258       >=10
17      25778.524       >=10
18      28344.260       >=10
19      28661.302       >=10
20      26028.344       >=10
21      20910.899       >=10
22      14863.931       >=10
23      9467.695        >=10
24      5502.058        >=10
25      2980.860        >=10
26      1544.432        >=10
27      1233.627        >=10
28      505.554 >=10
Groups below 10 expected recombinations: 0
```
