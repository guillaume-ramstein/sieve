cd /media/shd-sieve/N1/PEER
for n in 1 2 3 4 6 7 8 9 11 12 13 14 15; do
    ./bin/peertool -f peer.expression.csv -n $n -o peer_n-$n -i 1000
done