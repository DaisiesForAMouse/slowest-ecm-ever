import subprocess
import sys
import re
import csv
import time


def get_factors(num, cores):
    start = time.time()
    proc = subprocess.run(["./moduli", str(num), str(cores)],  capture_output = True)
    sol = proc.stdout.decode("utf-8").strip()
    print(sol)
    sol_re = re.search(r"SOLUTION: \d+ = (\d+) \* (\d+)", sol)
    p = sol_re.group(1);
    q = sol_re.group(2);
    return int(p), int(q), time.time() - start


def main():
    problems = []
    with open(sys.argv[1], newline="") as file:
            moduli = csv.reader(file, delimiter=",")
            for row in moduli:
                title = row[0]
                num = int(row[1])
                problems += [(title, num)]

    problems = sorted(problems, key = lambda p: p[1])
    with open(sys.argv[2], "w", newline="") as out_file:
        for title, num in problems:
            sol_writer = csv.writer(out_file, delimiter=",")
            p, q, time = get_factors(num, 8)
            time = round(time, 4)
            print("time spent factoring {}: {}s".format(num, time))
            sol_writer.writerow([title, num, p, q, time])
            out_file.flush()


if __name__ == "__main__":
    main()
    # get_factors(433859, 4);