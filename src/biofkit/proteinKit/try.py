a = 8

try:
    if (a == 8):
        raise (TypeError)
except TypeError as e:
    print('TypeError')
    raise

print('still alive')
