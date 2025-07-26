# import modules
from time import perf_counter
from random import randint
from sys import argv, exit, modules
import numpy as np
from sympy.ntheory import factorint
from memory_profiler import memory_usage
import importlib.util


def test(src, dst):
    # read data from csv file
    data = read_csv(src)
    results = data[:]

    for i, instance in enumerate(data):
        # test on instance
        time_data, runs = test_instance(instance)
        mean, sd = summarise_data(time_data)

        # add to results
        results[i] += [mean, sd, runs]
    # write to csv
    write_csv(results, dst)


def test_instance(instance, time_lim=100, min_runs=30, max_runs=500):
    # init. variables
    data, runs = [], 0
    p, g, n = instance[:3]
    fact = factorint(n)

    start = perf_counter()
    while ((perf_counter() - start < time_lim) and (runs < max_runs)) or (runs < min_runs):
        # random target
        h = pow(g, randint(0, n - 1), p)
        args = (p, g, h, n, fact) if pohlig else (p, g, h, n)

        # run function
        if metric:
            # test execution time
            s = perf_counter()
            x = f(*args)
            elapsed = perf_counter() - s
            data.append(elapsed)
        else:
            # test memory usage
            peak_mem, x = memory_usage((f, args), max_usage=True, retval=True)
            data.append(peak_mem)

        # verify correctness
        if x is None or h != pow(g, x, p):
            print(f'ERROR: incorrect for {p},{g},{h} - returned: {x}')
        runs += 1

    print(f'SUCCESS: performed {runs} runs on {p},{g}')
    return data, runs


def read_csv(src_path, delim=','):
    data = []
    # read from file
    with open(src_path, 'r') as f:
        for line in f:
            data.append([int(float(val)) if float(val).is_integer() else float(val) 
             for val in line.strip().split(delim)])

    # alert success and output data 
    print('SUCCESS: data read from file')
    return data
            

def write_csv(data, dst_path, delim=','):
    # convert data to string format
    output = ""
    for row in data:
        output += delim.join(map(str, row)) + '\n'
    
    # write to file
    with open(dst_path, 'w') as f:
        f.write(output)
    print('SUCCESS: data written to file')


def summarise_data(data):
    # mean of data, mu
    mu = np.mean(data)
    # sample standard deviation, s
    sigma = np.std(data, ddof=1)
    return mu, sigma


if __name__ == "__main__":
    if len(argv) != 8:
        print("Usage: python3 script.py <data_path> <results_path> <file_path> <module_name> <function_name> <metric> <testing_pohlig?>")
        exit()

    # parse arguments
    src, dst, file, module, func, metric, pohlig = argv[1:]
    metric = True if metric != "memory" and metric != "m" and metric != "False" and metric != "mem" else False
    pohlig = False if pohlig != "True" and pohlig != "1" else True

    # import dlp algorithm
    spec = importlib.util.spec_from_file_location(module, file)
    mod = importlib.util.module_from_spec(spec)
    modules[module] = mod
    spec.loader.exec_module(mod)
    f = getattr(mod, func)

    test(src, dst)
