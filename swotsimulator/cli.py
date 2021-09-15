import sys
import swotsimulator.mod_tools as mod_tools
import swotsimulator.settings as settings
import argparse
import logging


def run_swot_script():
    """Run SWOT Simulator"""
    import swotsimulator.run_simulator as run_simulator
    # Setup logging
    main_logger = logging.getLogger()
    main_logger.handlers = []
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    main_logger.addHandler(handler)
    main_logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument('params_file', nargs='?', type=str, default=None,
                        help='Path of the parameters file')
    parser.add_argument('--die-on-error', action='store_true', default=False,
                        help='Force simulation to quit on first error')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Display debug log messages')

    args = parser.parse_args()

    if args.params_file is None:
        logger.error('Please specify a parameter file')
        sys.exit(1)

    if args.debug is True:
        main_logger.setLevel(logging.DEBUG)

    file_param = args.params_file

    #p = mod_tools.load_python_file(file_param)
    p = settings.Parameters(settings.eval_config_file(file_param))
    try:
        run_simulator.run_simulator(p, args.die_on_error)
    except KeyboardInterrupt:
        logger.error('\nInterrupted by user (Ctrl+C)')
        sys.exit(1)
    sys.exit(0)


def run_nadir_script():
    """Run Altimeter Simulator"""
    import swotsimulator.run_simulator as run_simulator

    # Setup logging
    main_logger = logging.getLogger()
    main_logger.handlers = []
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    main_logger.addHandler(handler)
    main_logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument('params_file', nargs='?', type=str, default=None,
                        help='Path of the parameters file')
    parser.add_argument('--die-on-error', action='store_true', default=False,
                        help='Force simulation to quit on first error')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Display debug log messages')

    args = parser.parse_args()

    if args.params_file is None:
        logger.error('Please specify a parameter file')
        sys.exit(1)

    if args.debug is True:
        main_logger.setLevel(logging.DEBUG)

    file_param = args.params_file

    p = mod_tools.load_python_file(file_param)
    try:
        run_simulator.run_simulator(p, die_on_error=args.die_on_error,
                                    nadir_alone=True)
    except KeyboardInterrupt:
        logger.error('\nInterrupted by user (Ctrl+C)')
        sys.exit(1)
    sys.exit(0)
