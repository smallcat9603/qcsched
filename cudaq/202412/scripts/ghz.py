
import cudaq

@cudaq.kernel
def ghz(num_qubits: int):
    q = cudaq.qvector(num_qubits)
    h(q[0])
    for i in range(1, num_qubits):
        cx(q[i-1], q[i])

print(cudaq.draw(ghz, 10))

import numpy as np
print(np.array(cudaq.get_state(ghz, 10)))
