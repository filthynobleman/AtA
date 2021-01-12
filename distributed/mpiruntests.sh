mkdir tests


mkdir tests/10kx10k
for i in {1..8}
do
    printf "mpirun -np ${i} ./main 10 10000 10000 tests/10kx10k/times.${i}.txt\t"
    mpirun -np ${i} ./main 10 10000 10000 tests/10kx10k/times.${i}.txt
    printf "Done\n"
done


mkdir tests/20kx20k
for i in {1..8}
do
    printf "mpirun -np ${i} ./main 10 20000 20000 tests/20kx20k/times.${i}.txt\t"
    mpirun -np ${i} ./main 10 20000 20000 tests/20kx20k/times.${i}.txt
    printf "Done\n"
done


mkdir tests/10kx40k
for i in {1..8}
do
    printf "mpirun -np ${i} ./main 10 10000 40000 tests/10kx40k/times.${i}.txt\t"
    mpirun -np ${i} ./main 10 10000 40000 tests/10kx40k/times.${i}.txt
    printf "Done\n"
done


mkdir tests/40kx10k
for i in {1..8}
do
    printf "mpirun -np ${i} ./main 10 40000 10000 tests/40kx10k/times.${i}.txt\t"
    mpirun -np ${i} ./main 10 40000 10000 tests/40kx10k/times.${i}.txt
    printf "Done\n"
done


mkdir tests/5kx60k
for i in {1..8}
do
    printf "mpirun -np ${i} ./main 10 5000 60000 tests/5kx60k/times.${i}.txt\t"
    mpirun -np ${i} ./main 10 5000 60000 tests/5kx60k/times.${i}.txt
    printf "Done\n"
done


mkdir tests/60kx5k
for i in {1..8}
do
    printf "mpirun -np ${i} ./main 10 60000 5000 tests/60kx5k/times.${i}.txt\t"
    mpirun -np ${i} ./main 10 60000 5000 tests/60kx5k/times.${i}.txt
    printf "Done\n"
done