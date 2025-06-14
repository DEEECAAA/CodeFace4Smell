# This file is part of Codeface. Codeface is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, version 2.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Copyright 2013 by Siemens AG
# All Rights Reserved.

from tempfile import mkdtemp, NamedTemporaryFile
from shutil import rmtree
from collections import namedtuple
from subprocess import check_call
from textwrap import dedent
from os import getcwd, chdir, listdir, unlink, getenv, makedirs, environ
from os.path import split as pathsplit, join as pathjoin, isdir, exists, basename
from datetime import datetime
from time import strptime
from random import Random
import re
import yaml

_Author = namedtuple("_Author", ["name", "email"])
class Author(_Author):
    def __str__(self):
        return "{a.name} <{a.email}>".format(a=self)
Commit = namedtuple("Commit", ["author", "committer", "datetime", "filetree", "signoff", "tags"])
Tag = namedtuple("Tag", ["author", "datetime", "type"])
iso8601 = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\+\d{2}:\d{2})?")
iso8601_simple = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$")

class GitProject(object):
    '''
    Represents a Git repository with utility methods to create fake histories
    with fake authors, committers, rc and release tags.
    To get a physical git repository, use an instance of ExampleProject
    as a context manager using the with statement. This ensures that the git
    repository is properly deleted at the end of testing.
    '''
    def __init__(self, tagging="tag", randomise_email_case=False):
        '''Creates a repository with no commits'''
        self._authors = []
        self._commits = []
        self._tagging = tagging
        self._mboxes = {}
        self._randomise_email_case = randomise_email_case
        self._random = Random()
        self._random.seed(42)

    def __enter__(self):
        '''
        This function is called when entering a with statement.
        Here a real git repository is created and populated with the
        information in the project.
        '''
        self.directory = mkdtemp(prefix="codeface_test_project")
        cwd = getcwd()
        try:
            chdir(self.directory)

            def git(cmds, committer=None, commitdate=None, authordate=None):
                env = dict(environ)
                if committer:
                    env["GIT_COMMITTER_NAME"] = committer.name
                    env["GIT_COMMITTER_EMAIL"] = committer.email
                if commitdate:
                    env["GIT_COMMITTER_DATE"] = commitdate
                if authordate:
                    env["GIT_AUTHOR_DATE"] = authordate
                check_call(["git"] + cmds, env=env)

            git(["init"])
            next_release = 0
            next_rc = 0
            release_tags = []
            rc_tags = {}
            for i, c in enumerate(self._commits):
                for f in listdir("."):
                    if f != ".git":
                        if exists(f) and isdir(f):
                            rmtree(f)
                        elif exists(f):
                            unlink(f)
                for f, content in c.filetree.items():
                    dn, fn = pathsplit(f)
                    if dn and not exists(dn):
                        makedirs(dn)
                    with open(f, "w") as fd:
                        fd.write(content)
                git("add -A .".split())
                commitmsg = f"Commit {i}\n\nCommit message\n\n"
                for signer in c.signoff:
                    commitmsg += f"Signed-off-by: {str(signer)}\n"
                if not iso8601.match(c.datetime):
                    raise ValueError("expected iso8601 date (timezone is optional and set to +01:00 if not available)")
                timezoned_date = c.datetime + "+01:00" if iso8601_simple.match(c.datetime) else c.datetime
                git(["commit", "--author", str(c.author), "--date", timezoned_date, "-m", commitmsg],
                    committer=c.committer, commitdate=timezoned_date, authordate=timezoned_date)
                for tag in c.tags:
                    name = f"v{next_release}_{tag.type}"
                    if tag.type == "rc":
                        name += f"_{next_rc}"
                        rc_tags.setdefault(next_release, name)
                        next_rc += 1
                    elif tag.type == "release":
                        release_tags.append(name)
                        next_release += 1
                    git(["tag", name])

            # Applica i sottosistemi definiti
            if hasattr(self, "_subsystems"):
                subsys_path = pathjoin(self.directory, ".git", "subsystems.txt")
                with open(subsys_path, "w") as f:
                    for filepath, subsystem in self._subsystems.items():
                        f.write(f"{filepath}\t{subsystem}\n")

            # Allineamento delle due liste
            rcs = [rc_tags.get(i, release_tags[i]) for i in range(len(release_tags))]
            if len(release_tags) != len(rcs):
                raise ValueError(f"Mismatch between revisions ({len(release_tags)}) and rcs ({len(rcs)})")

            basename_dir = basename(self.directory)

            config_data = {
                "project": basename_dir,
                "repo": basename_dir,
                "description": f"{basename_dir} Description",
                "mailinglists": [
                    {"name": f"{basename_dir}.dev1", "type": "dev", "source": "generated"},
                    {"name": f"{basename_dir}.dev2", "type": "dev", "source": "generated"},
                    {"name": f"{basename_dir}.user1", "type": "user", "source": "generated"},
                    {"name": f"{basename_dir}.user2", "type": "user", "source": "generated"},
                ],
                "revisions": release_tags,
                "rcs": rcs,
                "tagging": self._tagging
            }

            with open(self.codeface_conf, "w") as fd:
                yaml.dump(config_data, fd, sort_keys=False)

            for ml_name, ml_file in self.mboxes:
                with open(ml_file, "w") as fd:
                    fd.write(self.mbox_contents(ml_name))

        finally:
            chdir(cwd)

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None:
            rmtree(self.directory)
        else:
            print("Left directory '{}' for inspection.".format(self.directory))
        del self.directory

    def add_author(self, name, email):
        author = Author(name, email)
        self._authors.append(author)
        return author

    def commit(self, author, committer, datetime, filetree, signoff):
        self._commits.append(Commit(author, committer, datetime, filetree, signoff, []))

    def tag_rc(self, author, datetime):
        self._commits[-1].tags.append(Tag(author, datetime, "rc"))

    def tag_release(self, author, datetime):
        self._commits[-1].tags.append(Tag(author, datetime, "release"))

    @property
    def gitdir(self):
        return self.directory

    @property
    def codeface_conf(self):
        return pathjoin(self.directory, ".git", "testproject.conf")

    @property
    def name(self):
        return basename(self.directory)

    @property
    def mboxes(self):
        project = basename(self.directory)
        return [(name, pathjoin(self.directory, ".git",
                 ".".join((project, name, "mbox"))))
                for name in ("dev1", "dev2", "user1", "user2")]

    @property
    def authors(self):
        return self._authors

    def email(self, mlist, sender, date, subject, content):
        cdate = datetime(*strptime(date, "%Y-%m-%dT%H:%M:%S")[:6]).strftime("%a, %d %b %Y %H:%M:%S")
        if self._randomise_email_case:
            subject = ''.join(self._random.choice((x.lower(), x.upper())) for x in subject)
            content = ''.join(self._random.choice((x.lower(), x.upper())) for x in content)
        self._mboxes.setdefault(mlist, []).append(dedent(
        """
        From MAILER-DAEMON Thu Jul 18 13:48:48 2013
        Path: example.com!not-for-mail
        From: {sender}
        Newsgroups: gmane.codeface.test.project
        Subject: {subject}
        Date: {date}
        Approved: auto
        Message-ID: <{messageid}@example.com>
        NNTP-Posting-Host: example.com
        Mime-Version: 1.0
        Content-Type: text/plain; charset=us-ascii; format=flowed
        Content-Transfer-Encoding: 7bit
        X-Complaints-To: complaints@example.com
        NNTP-Posting-Date: {date}
        User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:0.9.8) Gecko/20020205
        X-Accept-Language: en-us
        Original-To: codeface.test.project@example.com
        Precedence: bulk
        X-Mailing-List: codeface.test.project@example.com

        {content}""").
        format(date=cdate, sender=sender, subject=subject,
               content=content, messageid=len(self._mboxes.get(mlist, []))))

    def mbox_contents(self, mlist):
        return ("\n\n".join(self._mboxes.get(mlist, [])) + "\n\n").lstrip()

    def assign_subsystem(self, filepath, subsystem):
        if not hasattr(self, "_subsystems"):
            self._subsystems = {}
        self._subsystems[filepath] = subsystem