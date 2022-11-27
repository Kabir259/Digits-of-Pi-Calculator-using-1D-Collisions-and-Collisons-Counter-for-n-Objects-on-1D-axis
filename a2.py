class heapNode:
    def __init__(self, time, ind, position, tlu):
        self.time = time
        self.ind = ind
        self.position = position
        self.tlu = tlu
        self.output = [self.time, self.ind, self.position]

class Heap:
    def __init__(self, m):
        self.a = []
        self.indexTracker = [-1] * m

    def size(self):
        return len(self.a)

    def isEmpty(self):
        return self.size() == 0

    def heapifyUp(self, index):
        childIndex = index
        while childIndex > 0:
            parentIndex = (childIndex - 1) // 2
            if self.a[childIndex].time < self.a[parentIndex].time:
                self.a[childIndex], self.a[parentIndex] = self.a[parentIndex], self.a[childIndex]
                self.indexTracker[self.a[childIndex].ind], self.indexTracker[self.a[parentIndex].ind] = \
                    self.indexTracker[self.a[parentIndex].ind], self.indexTracker[self.a[childIndex].ind]
                childIndex = parentIndex

            if self.a[childIndex].time == self.a[parentIndex].time:
                if self.a[childIndex].ind < self.a[parentIndex].ind:
                    self.a[childIndex], self.a[parentIndex] = self.a[parentIndex], self.a[childIndex]
                    self.indexTracker[self.a[childIndex].ind], self.indexTracker[self.a[parentIndex].ind] = \
                        self.indexTracker[self.a[parentIndex].ind], self.indexTracker[self.a[childIndex].ind]
                    childIndex = parentIndex
                else:
                    break
            else:
                break

    def insert(self, element):
        # element = heapNode(time, ind, position, tlu)
        self.a.append(element)
        self.indexTracker[element.ind] = self.size() - 1
        self.heapifyUp(self.size() - 1)

    def heapifyDown(self, index):
        parentIndex = index
        leftChildIndex = 2 * parentIndex + 1
        rightChildIndex = 2 * parentIndex + 2

        while leftChildIndex < self.size():
            lowestIndex = parentIndex

            if self.a[lowestIndex].time > self.a[leftChildIndex].time:
                lowestIndex = leftChildIndex
            if rightChildIndex < self.size() and self.a[lowestIndex].time > self.a[rightChildIndex].time:
                lowestIndex = rightChildIndex
            if lowestIndex == parentIndex:
                if self.a[lowestIndex].ind > self.a[leftChildIndex].ind:
                    lowestIndex = leftChildIndex
                if rightChildIndex < self.size() and self.a[lowestIndex].ind > self.a[rightChildIndex].ind:
                    lowestIndex = rightChildIndex
                else:  # if lowestIndex == parentIndex even after comparing indices of both lowest with left/right child
                    break

            self.a[parentIndex], self.a[lowestIndex] = self.a[lowestIndex], self.a[parentIndex]
            self.indexTracker[self.a[lowestIndex].ind], self.indexTracker[self.a[parentIndex].ind] = \
                self.indexTracker[self.a[parentIndex].ind], self.indexTracker[self.a[lowestIndex].ind]

            parentIndex = lowestIndex
            leftChildIndex = 2 * parentIndex + 1
            rightChildIndex = 2 * parentIndex + 2

    def removeMin(self):
        first = self.a[0]
        last = self.a[-1]
        self.a[0] = heapNode(last.time,last.ind,last.position,0)
        self.indexTracker[first.ind]=-1
        popped = self.a.pop()
        self.indexTracker[popped.ind] = 0
        self.heapifyDown(0)
        return first.output

    def top(self):
        return self.a[0]

    def modify(self, Node):
        if self.indexTracker[Node.ind] == -1:
            return
        childIndex = self.indexTracker[Node.ind]
        parentIndex = (childIndex - 1) // 2
        self.a[childIndex] = Node
        self.heapifyDown(childIndex) if self.a[parentIndex].time > self.a[childIndex].time else self.heapifyUp(childIndex)


def simpleColltime(x1,x2,v1,v2,T):
    return T + 1 if v2 >= v1 else abs((x2 - x1)/(v2 - v1))


def collVel(u1, u2, m1, m2):
    if m1 == float('inf'):
        v1 = u1
        v2 = -u2
        return [v1, v2]
    v1 = ((m1 - m2) / (m1 + m2)) * u1 + 2 * m2 * u2 / (m1 + m2)
    v2 = ((m2 - m1) / (m1 + m2)) * u2 + 2 * m1 * u1 / (m1 + m2)
    return [v1, v2]


def newNodeMaker(x1, x2, v2, v1, tlu, t0, T):
    if v2 >= v1:
        time = T + 1
        position = x2 + v2 * time
        return [time, position]
    time = abs((x2 - x1 - v1 * (tlu - t0)) / (v2 - v1)) + tlu
    position = x2 + v2 * time
    return [time, position]


def listCollisions(M, x, v, m, T):
    finalOutput = []
    numCollisions = len(finalOutput)
    netTime = 0
    n = len(M)
    heap = Heap(n)

    for i in range(n - 1):
        time = simpleColltime(x[i],x[i+1],v[i],v[i+1],T)
        heap.insert(heapNode(time, i, x[i] + v[i] * time, 0))

    placeholder = [0] * n

    while numCollisions < m and netTime < T:
        trueOutput = []
        track = []
        rootTime = heap.top().time
        tlu = rootTime

        if tlu >= T:
            break
        else:
            while not heap.size() == 0 and numCollisions <= m:
                if heap.top().time == rootTime:
                    extractedRoot = heap.removeMin()
                    track.extend([extractedRoot[1]])
                    trueOutput.extend([(extractedRoot[0], extractedRoot[1], extractedRoot[2])])
                    finalOutput.extend([(round(extractedRoot[0], 4), extractedRoot[1], round(extractedRoot[2], 4))])
                    numCollisions += 1
                else:
                    break

        t = len(track)

        for i in range(t):
            v[track[i]] = ((M[track[i]] - M[track[i] + 1]) * v[track[i]] + 2 * M[track[i] + 1] * v[track[i] + 1]) / (M[track[i]] + M[track[i] + 1])
            v[track[i] + 1] = ((M[track[i] + 1] - M[track[i]]) * v[track[i] + 1] + 2 * M[track[i]] * v[track[i]]) / (M[track[i]] + M[track[i] + 1])
            placeholder[track[i]] = placeholder[track[i] + 1] = tlu
            x[track[i]] = x[track[i] + 1] = trueOutput[i][2]

            if track[i] > 0:
                time1 = abs(x[track[i]] - x[track[i] - 1] - v[track[i] - 1] * (tlu - placeholder[track[i] - 1])) / abs(v[track[i]] - v[track[i] - 1])
                position = x[track[i]] + v[track[i]] * time1
                heap.modify(heapNode(T + tlu, track[i] - 1, x[track[i]] + v[track[i]] * T, 0)) if v[track[i]] >= v[track[i] - 1] else heap.modify(heapNode(time1 + tlu, track[i] - 1, position, 0))

            if track[i] < len(M) - 2:
                if v[track[i] + 2] >= v[track[i] + 1]:
                    time2 = T
                    position = x[track[i] + 1] + v[track[i] + 1] * time2
                    heap.modify(heapNode(T + tlu, track[i] + 1, position, 0))
                else:
                    time2 = abs(x[track[i] + 2] + v[track[i] + 2] * (tlu - placeholder[track[i] + 2]) - x[track[i] + 1]) / abs(v[track[i] + 2] - v[track[i] + 1])
                    position = x[track[i] + 1] + v[track[i] + 1] * time2
                    heap.modify(heapNode(time2 + tlu, track[i] + 1, position, 0))


        for i in range(len(trueOutput)):
            newnode = heapNode(T, trueOutput[i][1], 0, 0)
            heap.insert(newnode)
        if numCollisions >= m:
            break
    return finalOutput


print(listCollisions([10000.0, 1.0, 100.0], [0.0, 1.0, 2.0], [0.0, 0.0, -1.0], 6, 10.0))
