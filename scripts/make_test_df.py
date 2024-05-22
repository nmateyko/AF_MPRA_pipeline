import pandas as pd
import random

data = []

for i in range(1000000):
    data.append({
        'Read ID': f'this_is_sequence_{i}',
        'Sequence': "".join(random.choices("ACGT", k=105)),
        'Quality': ''.join(random.choices("ABCDEF", k=105)),
        'CIGAR': ''.join(random.choices("0123456789MIDNSHP=X", k=10)),
        'Flag': 'flag',
        'Reference ID': 'reference_id',
        'Position': 234,
        'Probe_name': 'probe_name',
        'barcode': "".join(random.choices("ACGT", k=20))
    })

df = pd.DataFrame(data)
print(df.memory_usage(deep=True))