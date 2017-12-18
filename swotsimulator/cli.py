import sys
import swotsimulator.mod_tools as mod_tools
import argparse
import logging
logger = logging.getLogger()
handler = logging.StreamHandler()


def run_swot_script():
    """Run SWOT Simulator"""
    import swotsimulator.run_simulator as run_simulator
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Display debug log messages')
    args = parser.parse_args()
    if args.debug is True:
        logger.setLevel(logging.DEBUG)
    '''
    if len(sys.argv) < 2:
        logger.error('Please specify a parameter file')
        sys.exit(1)
    else:
        file_param = str(sys.argv[1])

    p = mod_tools.load_python_file(file_param)
    run_simulator.run_simulator(p)

def run_nadir_script():
    """Run Altimeter Simulator"""
    import swotsimulator.run_simulator as run_simulator
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    if len(sys.argv) < 2:
        logger.error('Please specify a parameter file')
        sys.exit(1)
    else:
        file_param = str(sys.argv[1])

    p = mod_tools.load_python_file(file_param)
    run_simulator.run_nadir(p)
    
