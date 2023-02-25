import os
import shutil
def main():
    # remove all folders whose name starts with "temp"
    for folder in os.listdir():
        if folder.startswith("temp"):
            shutil.rmtree(folder)

if __name__ == "__main__":
    main()