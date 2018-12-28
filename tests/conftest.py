from snakescale import WRAPPER_ROOT_PATH, DEFAULT_SCALE_FILENAME

collect_ignore = ['setup.py']

print(WRAPPER_ROOT_PATH)

for path in WRAPPER_ROOT_PATH.rglob(DEFAULT_SCALE_FILENAME):
    collect_ignore.append(str(path))

print(collect_ignore)