import sys

# https://stackoverflow.com/a/51704613
try:
    from pip import main as pipmain
except ImportError:
    from pip._internal import main as pipmain

if __name__ == '__main__':
    try:
        pipmain(['install', '-r', 'requirements.txt', '--upgrade', '--upgrade-strategy=eager'])
    except Exception as e:
        print("Encountered exception: ", e)
        print("Press enter to close.")
        input()