import threading
import time

def worker(event):
    t = 0
    while True:
        t += 1
        print(f'{t=}')
        time.sleep(1)
        event.wait()
        if t == 10:
            break
        
event = threading.Event()
thread = threading.Thread(target=worker, args=(event,))
thread.start()

event.set()

time.sleep(3)
print('stop thread')
event.clear()

time.sleep(3)
print('resume thread')
event.set()