import os
import os.path
import sqlite3

# removes the database at the start
if os.path.exists("courses.db"):
    os.remove("courses.db")

db = sqlite3.connect("courses.db")
db.isolation_level = None

# creates required tables for the database
def create_tables():
    db.executescript("""
        CREATE TABLE Teachers
        (
            id              INTEGER PRIMARY KEY,
            name            TEXT
        );
        CREATE TABLE Students
        (
            id              INTEGER PRIMARY KEY,
            name            TEXT
        );
        CREATE TABLE Courses
        (
            id              INTEGER PRIMARY KEY,
            name            TEXT,
            credits         INTEGER
        );
        CREATE TABLE CourseTeachers
        (
            id             INTEGER PRIMARY KEY,
            course_id      TEXT REFERENCES Courses(id),
            teacher_id     TEXT REFERENCES Teachers(id)
        );
        CREATE TABLE StudentCredits
        (
            id             INTEGER PRIMARY KEY,
            student_id     TEXT REFERENCES Students(id),
            course_id      TEXT REFERENCES Courses(id),
            date           DATE,
            grade          INTEGER
        );
        CREATE TABLE Groups
        (
            id             INTEGER PRIMARY KEY,
            name           TEXT
        );
        CREATE TABLE GroupTeachers
        (
            id             INTEGER PRIMARY KEY,
            group_id       INTEGER REFERENCES Groups(id),
            teacher_id     TEXT REFERENCES Teachers(id)
        );
        CREATE TABLE GroupStudents
        (
            id             INTEGER PRIMARY KEY,
            group_id       INTEGER REFERENCES Groups(id),
            student_id     TEXT REFERENCES Students(id)
        );          
        """
    )

# adds a teacher to the database
def create_teacher(name):
    sql = """INSERT INTO Teachers(name) VALUES(?)"""
    id = db.execute(sql, [name]).lastrowid
    return id


# adds a course to the database
def create_course(name, credits, teacher_ids):
    sql = """INSERT INTO Courses(name, credits) VALUES(?,?)"""
    id = db.execute(sql, [name, credits]).lastrowid
    sql = """INSERT INTO CourseTeachers(course_id, teacher_id) VALUES(?,?)"""
    db.executemany(sql, [(id,t_id) for t_id in teacher_ids])
    return id


# adds a student to the database
def create_student(name):
    sql = """   INSERT INTO Students(name) VALUES(?)"""
    id = db.execute(sql, [name]).lastrowid
    return id


# adds credits to a student
def add_credits(student_id, course_id, date, grade):
    sql = """   INSERT INTO StudentCredits(student_id, course_id, date, grade) VALUES(?,?,?,?)"""
    id = db.execute(sql, [student_id, course_id, date, grade]).lastrowid
    return id


# adds a group to the database
def create_group(name, teacher_ids, student_ids):
    sql = """INSERT INTO Groups(name) VALUES(?)"""
    id = db.execute(sql, [name]).lastrowid
    sql = """INSERT INTO GroupTeachers(group_id, teacher_id) VALUES(?,?)"""
    db.executemany(sql, [(id,t_id) for t_id in teacher_ids])
    sql = """INSERT INTO GroupStudents(group_id, student_id) VALUES(?,?)"""
    db.executemany(sql, [(id,s_id) for s_id in student_ids])
    return id


# looks for courses where a specific teacher is teaching (in alphabetical order)
def courses_by_teacher(teacher_name):
    sql = """
        SELECT C.name
        FROM CourseTeachers AS CT
        JOIN Courses AS C
            ON CT.course_id = C.id
        JOIN Teachers AS T
            ON CT.teacher_id = T.id
        WHERE T.name = ?
        ORDER BY C.name
        """
    res = db.execute(sql, [teacher_name]).fetchall()
    return [r[0] for r in res]


# gets the amount of credits a teacher has given
def credits_by_teacher(teacher_name):
    sql = """
        SELECT SUM(C.credits)
        FROM Teachers AS T, CourseTeachers AS CT, Courses AS C, StudentCredits AS SC
        WHERE T.id = CT.teacher_id
            AND C.id = CT.course_id
            AND C.id = SC.course_id
            AND T.name = ?
        """
    res = db.execute(sql, [teacher_name]).fetchone()
    return res[0]



# gets courses which a student has completed with grades included
def courses_by_student(student_name):
    sql = """
        SELECT C.name, SC.grade
        FROM Students AS S, Courses AS C, StudentCredits AS SC
        WHERE S.id = SC.student_id
            AND C.id = SC.course_id
            AND S.name = ?
        ORDER BY C.name
        """
    res = db.execute(sql, [student_name]).fetchall()
    return res


# gets the amount of credits from a specific year
def credits_by_year(year):
    sql = """
        SELECT SUM(credits)
        FROM StudentCredits AS SC, Courses AS C
        WHERE C.id = SC.course_id
            AND SUBSTR(date,1,4) = ?
        """
    res = db.execute(sql, [str(year)]).fetchone()
    return res[0]


# gets a grade distribution for a course (grades in order 1-5)
def grade_distribution(course_name):
    tulos = {n:0 for n in range(1,6)}
    sql = """
        SELECT grade, COUNT(*)
        FROM Courses AS C
        LEFT JOIN StudentCredits AS SC
        WHERE C.id = SC.course_id
            AND C.name = ?
        GROUP BY grade
        """
    res = db.execute(sql, [course_name]).fetchall()
    return tulos | dict(res)


# gets list of courses (name, # of teachers, # of students) in alphabetical order
def course_list():
    sql = """
    SELECT C.name,
            (SELECT IFNULL(COUNT(A.course_id),0)
            FROM CourseTeachers AS A WHERE C.id = A.course_id),
            (SELECT IFNULL(COUNT(B.course_id),0)
            FROM StudentCredits AS B WHERE C.id = B.course_id)
    FROM Courses AS C
    ORDER BY C.name
    """
    res = db.execute(sql).fetchall()
    return res


# gets a list of teachers and their courses (teachers and courses in alphabetical order)
def teacher_list():
    sql = """
    SELECT T.name, GROUP_CONCAT(C.name)
    FROM Teachers AS T, Courses AS C, CourseTeachers AS CT
    WHERE T.id = CT.teacher_id
        AND C.id = CT.course_id
    GROUP BY T.name
    ORDER BY T.name, C.name
    """
    res = db.execute(sql).fetchall()

    res = [(r[0], r[1].split(",")) for r in res]
    return res


# gets all the students in a group (in alphabetical order)
def group_people(group_name):
    sql = """
    SELECT T.name
    FROM Groups AS G, Teachers AS T, GroupTeachers AS GT
    WHERE G.id = GT.group_id
        AND T.id = GT.teacher_id
        AND G.name = ?
    UNION
    SELECT S.name
    FROM Groups AS G, Students AS S, GroupStudents AS GS
    WHERE G.id = GS.group_id
        AND S.id = GS.student_id
        AND G.name = ?
    """
    res = db.execute(sql, [group_name, group_name]).fetchall()

    res = [r[0] for r in res]
    return res


# gets the # of credits for each group (in alphabetical order)
def credits_in_groups():
    sql = """
    SELECT G.name, IFNULL(SUM(C.credits),0)
    FROM Groups AS G
    LEFT JOIN GroupStudents AS GS
        ON G.id = GS.group_id
    LEFT JOIN StudentCredits AS SC
        ON SC.student_id = GS.student_id
    LEFT JOIN Courses AS C
        ON C.id = SC.course_id
    GROUP BY G.name
    ORDER BY G.name
    """
    res = db.execute(sql).fetchall()

    return res


# gets the groups which has a specific teacher and a student (in alphabetical order)
def common_groups(teacher_name, student_name):
    sql = """
    SELECT G.name
    FROM Groups AS G, Teachers AS T, Students AS S
    LEFT JOIN GroupStudents AS GS
        ON G.id = GS.group_id
    LEFT JOIN GroupTeachers AS GT
        ON G.id = GT.group_id
    WHERE T.id = GT.teacher_id
        AND S.id = GS.student_id
        AND T.name = ?
        AND S.name = ?
    ORDER BY G.name
    """
    res = db.execute(sql, [teacher_name, student_name]).fetchall()

    res = [r[0] for r in res]
    return res