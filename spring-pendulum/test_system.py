import system as s


def main():
    num = int(input("Number of pendulums: "))
    T = int(input("Duration (in seconds): "))
    print("Danke!")
    
    a = s.System(num, 1, 9.8, 0.3, 1, 0.05, T)
    a.animate()
    return None
    
    
if __name__ == '__main__':
    main()
