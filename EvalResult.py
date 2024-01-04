from PatrolRobot import search
from argparse import ArgumentParser
from collections import Counter
import yaml
from tqdm import tqdm

def eval(args):
    with open(args.param) as file:
        param = yaml.safe_load(file)

    success = 0
    pathLen = 0
    for _ in tqdm(range(args.num_trials)):
        patrols, catchers = search(param, verbose=False)
        if len(catchers) > 0:
            epochPathLen = 0
            for name in catchers:
                epochPathLen += patrols[name].pathLen
            epochPathLen /= len(catchers)
            pathLen += epochPathLen
        numCatchs = Counter([p.success for p in patrols.values()])[True]
        if numCatchs == args.num_targets:
            success += 1
    print(f"detectionRate {100 * success/args.num_trials :.3f}, avgPathLen {pathLen/args.num_trials}:.3f")



if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--param', type=str, default='param.yaml')
    parser.add_argument('--num-targets', type=int, default=2)
    parser.add_argument('--num-trials', type=int, default=10)

    args = parser.parse_args()
    eval(args)
