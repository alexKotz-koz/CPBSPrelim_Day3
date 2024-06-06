import pandas as pd
import numpy as np
import logging
import argparse
import time
import os

from components.importSample import importsample

logDir = "data/logs"
os.makedirs(logDir, exist_ok=True)

logging.basicConfig(
    filename=f"{logDir}/app.log",
    filemode="w",
    format="%(message)s",
    level=logging.INFO,
)


def main():
    print("Hello World!")
    importsample()


if __name__ == "__main__":
    main()
