import sys


def main(metrics_file, output_file):
    with open(metrics_file) as f:
        while not f.readline().startswith('MEDIAN_INSERT_SIZE'):
            continue
        line = f.readline()
        line_split = line.split('\t')
        mean_insert_size = line_split[4]
    with open(output_file, 'w') as f:
        f.write(mean_insert_size)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
